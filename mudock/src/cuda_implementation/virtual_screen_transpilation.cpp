#include <alloca.h>
#include <cstddef>
#include <cstring>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/cuda_implementation/virtual_screen.cuh>
#include <mudock/grid.hpp>
#include <mudock/utils.hpp>
#include <span>

// Keep it to 32 to enable warp optimizations
#define BLOCK_SIZE 32

namespace mudock {
  static constexpr std::size_t max_non_bonds{1024 * 10};

  void init_texture_memory(const grid_map& map, cuda_object<fp_type>& tex_obj) {
    tex_obj.alloc(map.index.size_x() * map.index.size_y() * map.index.size_z());
    tex_obj.copy_host2device(map.data());
  }

  void set_device(const std::size_t gpu_id);

  void call_kernel(const int batch_ligands,
                   const int num_generations,
                   const int tournament_length,
                   const fp_type mutation_prob,
                   const int chromosome_number,
                   const int chromosome_stride,
                   const int atom_stride,
                   const int rotamers_stride,
                   const int nonbond_stride,
                   const fp_type* __restrict__ original_ligand_x,
                   const fp_type* __restrict__ original_ligand_y,
                   const fp_type* __restrict__ original_ligand_z,
                   fp_type* __restrict__ scratch_ligand_x,
                   fp_type* __restrict__ scratch_ligand_y,
                   fp_type* __restrict__ scratch_ligand_z,
                   const fp_type* __restrict__ ligand_vol,
                   const fp_type* __restrict__ ligand_solpar,
                   const fp_type* __restrict__ ligand_charge,
                   const int* __restrict__ ligand_num_hbond,
                   const fp_type* __restrict__ ligand_Rij_hb,
                   const fp_type* __restrict__ ligand_Rii,
                   const fp_type* __restrict__ ligand_epsij_hb,
                   const fp_type* __restrict__ ligand_epsii,
                   const int* __restrict__ ligand_num_nonbonds,
                   const int* __restrict__ ligand_nonbond_a1,
                   const int* __restrict__ ligand_nonbond_a2,
                   const int* __restrict__ ligand_num_atoms,
                   const int* __restrict__ ligand_num_rotamers,
                   const int* __restrict__ ligand_fragments,
                   const int* __restrict__ frag_start_atom_index,
                   const int* __restrict__ frag_stop_atom_index,
                   chromosome* __restrict__ chromosomes,
                   const fp_type minimum_x,
                   const fp_type minimum_y,
                   const fp_type minimum_z,
                   const fp_type maximum_x,
                   const fp_type maximum_y,
                   const fp_type maximum_z,
                   const fp_type center_x,
                   const fp_type center_y,
                   const fp_type center_z,
                   const int index_n_x,
                   const int index_n_xy,
                   const fp_type* const __restrict__* const __restrict__ atom_textures,
                   const int* __restrict__ atom_tex_indexes,
                   const fp_type* electro_texture,
                   const fp_type* desolv_texture,
                   XORWOWState* __restrict__ state,
                   fp_type* __restrict__ ligand_scores,
                   chromosome* __restrict__ best_chromosomes);

  virtual_screen_cuda::virtual_screen_cuda(const knobs k,
                                                     const std::size_t gpu_id,
                                                     std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                                                     std::shared_ptr<const grid_map>& electro_map,
                                                     std::shared_ptr<const grid_map>& desolv_map)
      : configuration(k),
        center(electro_map.get()->center),
        minimum_coord(electro_map.get()->minimum_coord),
        maximum_coord(electro_map.get()->maximum_coord),
        index_map(electro_map.get()->index) {
    // Set device
    set_device(gpu_id);
    // TODO move this part into the cuda_worker -> once per GPU
    // Allocate grid maps
    assert(electro_map.get()->index == desolv_map.get()->index &&
           desolv_map.get()->index == grid_atom_maps.get()->get_index());
    init_texture_memory(*electro_map.get(), electro_tex);

    init_texture_memory(*desolv_map.get(), desolv_tex);

    atom_texs.wrappers_pointer.alloc(num_device_map_textures());
    atom_texs.wrappers.reserve(num_device_map_textures());
    for (int index{0}; index < num_device_map_textures(); ++index) {
      atom_texs.wrappers.emplace_back();
      cuda_object<fp_type>& atom_tex = atom_texs.wrappers.back();
      const grid_map& grid_atom =
          grid_atom_maps.get()->get_atom_map(autodock_type_from_map(static_cast<device_map_textures>(index)));
      init_texture_memory(grid_atom, atom_tex);
      atom_texs.wrappers_pointer.host[index] = atom_tex.dev_pointer();
    }
    atom_texs.wrappers_pointer.copy_host2device();

    // Grid spacing fixed to 0.5 Angstrom
    // setup_constant_memory(electro_map.get()->minimum_coord,
    //                       electro_map.get()->maximum_coord,
    //                       electro_map.get()->center,
    //                       fp_type{2});
  }

  void virtual_screen_cuda::operator()(batch& incoming_batch) {
    const std::size_t batch_atoms    = incoming_batch.batch_max_atoms;
    const std::size_t batch_rotamers = incoming_batch.batch_max_rotamers;
    const std::size_t batch_ligands  = incoming_batch.num_ligands;
    // Resize data structures
    const std::size_t tot_atoms_in_batch = batch_ligands * batch_atoms;
    // const std::size_t tot_atoms_in_population     = tot_atoms_in_batch * configuration.population_number;
    const std::size_t tot_rotamers_atoms_in_batch = tot_atoms_in_batch * batch_rotamers;
    const std::size_t tot_rotamers_in_batch       = batch_ligands * batch_rotamers;
    const std::size_t batch_nonbonds              = batch_ligands * max_non_bonds;
    // Use double buffering on the GPU for actual and next population at each iteration
    const std::size_t population_stride = configuration.population_number * 2;
    original_ligand_x.alloc(tot_atoms_in_batch);
    original_ligand_y.alloc(tot_atoms_in_batch);
    original_ligand_z.alloc(tot_atoms_in_batch);
    scratch_ligand_x.alloc(tot_atoms_in_batch);
    scratch_ligand_y.alloc(tot_atoms_in_batch);
    scratch_ligand_z.alloc(tot_atoms_in_batch);
    ligand_fragments.alloc(tot_rotamers_atoms_in_batch);
    frag_start_atom_indices.alloc(tot_rotamers_in_batch);
    frag_stop_atom_indices.alloc(tot_rotamers_in_batch);
    ligand_vol.alloc(tot_atoms_in_batch);
    ligand_solpar.alloc(tot_atoms_in_batch);
    ligand_charge.alloc(tot_atoms_in_batch);
    ligand_Rij_hb.alloc(tot_atoms_in_batch);
    ligand_Rii.alloc(tot_atoms_in_batch);
    ligand_epsij_hb.alloc(tot_atoms_in_batch);
    ligand_epsii.alloc(tot_atoms_in_batch);
    // TODO check if it is required
    ligand_num_hbond.alloc(tot_atoms_in_batch);
    ligand_num_atoms.alloc(batch_ligands);
    ligand_num_rotamers.alloc(batch_ligands);
    ligand_scores.alloc(batch_ligands);
    best_chromosomes.alloc(batch_ligands);
    // Bonds
    num_nonbonds.alloc(batch_ligands);
    nonbond_a1.alloc(batch_nonbonds);
    nonbond_a2.alloc(batch_nonbonds);
    // GA Data structures
    chromosomes.alloc(population_stride * batch_ligands);
    // Support data precomputation
    map_texture_index.alloc(tot_atoms_in_batch);

    // Copy data
    std::size_t index{0};
    // TODO pragma OMP improve performance
    // Keep the fragments for the output
    std::vector<fragments<static_containers>> batch_fragments;
    batch_fragments.reserve(batch_ligands);
    // Support data structures
    grid<uint_fast8_t, index2D> nbmatrix{{static_cast<int>(batch_atoms), static_cast<int>(batch_atoms)}};
    std::vector<non_bond_parameter> non_bond_list;
    non_bond_list.reserve(max_non_bonds);
    for (auto& ligand: std::span(incoming_batch.molecules.data(), batch_ligands)) {
      const int stride_atoms = index * batch_atoms;
      // Atoms and bonds
      const int num_atoms                    = ligand.get()->num_atoms();
      ligand_num_atoms.host_pointer()[index] = num_atoms;
      // TODO bonds
      // Place the molecule to the center of the target protein
      const auto x = ligand.get()->get_x(), y = ligand.get()->get_y(), z = ligand.get()->get_z();
      const auto ligand_center_of_mass = compute_center_of_mass(x, y, z);
      translate_molecule(x,
                         y,
                         z,
                         center.x - ligand_center_of_mass.x,
                         center.y - ligand_center_of_mass.y,
                         center.z - ligand_center_of_mass.z);

      std::memcpy((void*) (original_ligand_x.host_pointer() + stride_atoms),
                  x.data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void*) (original_ligand_y.host_pointer() + stride_atoms),
                  y.data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void*) (original_ligand_z.host_pointer() + stride_atoms),
                  z.data(),
                  num_atoms * sizeof(fp_type));
      // Fragments
      // Find out the rotatable bonds in the ligand
      auto graph = make_graph(ligand.get()->get_bonds());
      // TODO check this assignment
      batch_fragments[index]  = {graph, ligand.get()->get_bonds(), ligand.get()->num_atoms()};
      const auto& l_fragments = batch_fragments[index];

      // Randomly initialize the population
      const auto num_rotamers                   = l_fragments.get_num_rotatable_bonds();
      ligand_num_rotamers.host_pointer()[index] = num_rotamers;
      const int stride_masks                    = index * batch_rotamers * batch_atoms;
      const int stride_rotamers                 = index * batch_rotamers;
      assert(batch_rotamers > ligand.get()->num_rotamers());
      for (int rot = 0; rot < num_rotamers; ++rot) {
        std::memcpy((void*) (ligand_fragments.host_pointer() + stride_masks + rot * batch_atoms),
                    l_fragments.get_mask(rot).data(),
                    num_atoms * sizeof(int));
        const auto [start_index, stop_index]                          = l_fragments.get_rotatable_atoms(rot);
        frag_start_atom_indices.host_pointer()[stride_rotamers + rot] = start_index;
        frag_stop_atom_indices.host_pointer()[stride_rotamers + rot]  = stop_index;
      }

      // Weed bonds
      // TODO can be accelerated on the GPU
      // No need to move data -> already computed on the GPU
      nbmatrix.reset();
      nonbonds(nbmatrix, ligand.get()->get_bonds(), num_atoms);
      non_bond_list.clear();
      weed_bonds(nbmatrix, non_bond_list, num_atoms, l_fragments);
      // if constexpr (is_debug())
      //   if (non_bond_list.size() >= max_non_bonds) {
      //     throw std::runtime_error("Bond list size exceed maximum value " +
      //                              std::to_string(non_bond_list.size()) + ".");
      //   }

      num_nonbonds.host_pointer()[index] = non_bond_list.size();
      const int stride_nonbonds          = index * max_non_bonds;
      int nonbond_index{0};
      for (auto& bond: non_bond_list) {
        nonbond_a1.host_pointer()[stride_nonbonds + nonbond_index] = bond.a1;
        nonbond_a2.host_pointer()[stride_nonbonds + nonbond_index] = bond.a2;
        ++nonbond_index;
      }

      // Autodock typing
      std::memcpy((void*) (ligand_vol.host_pointer() + stride_atoms),
                  ligand.get()->get_vol().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void*) (ligand_solpar.host_pointer() + stride_atoms),
                  ligand.get()->get_solpar().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void*) (ligand_charge.host_pointer() + stride_atoms),
                  ligand.get()->get_charge().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void*) (ligand_Rij_hb.host_pointer() + stride_atoms),
                  ligand.get()->get_Rij_hb().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void*) (ligand_Rii.host_pointer() + stride_atoms),
                  ligand.get()->get_Rii().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void*) (ligand_epsij_hb.host_pointer() + stride_atoms),
                  ligand.get()->get_epsij_hb().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void*) (ligand_epsii.host_pointer() + stride_atoms),
                  ligand.get()->get_epsii().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void*) (ligand_num_hbond.host_pointer() + stride_atoms),
                  ligand.get()->get_num_hbond().data(),
                  num_atoms * sizeof(int));

      std::size_t atom_index{0};
      for (auto& autodock_t: ligand.get()->get_autodock_type()) {
        const auto map_index = static_cast<int>(map_from_autodock_type(autodock_t));
        map_texture_index.host_pointer()[stride_atoms + atom_index] = map_index;
        ++atom_index;
      }

      ++index;
    }
    // Copy in
    ligand_num_atoms.copy_host2device();
    original_ligand_x.copy_host2device();
    original_ligand_y.copy_host2device();
    original_ligand_z.copy_host2device();
    ligand_fragments.copy_host2device();
    ligand_num_rotamers.copy_host2device();
    frag_start_atom_indices.copy_host2device();
    frag_stop_atom_indices.copy_host2device();
    ligand_vol.copy_host2device();
    ligand_solpar.copy_host2device();
    ligand_charge.copy_host2device();
    ligand_Rij_hb.copy_host2device();
    ligand_Rii.copy_host2device();
    ligand_epsij_hb.copy_host2device();
    ligand_epsii.copy_host2device();
    ligand_num_hbond.copy_host2device();
    map_texture_index.copy_host2device();
    num_nonbonds.copy_host2device();
    nonbond_a1.copy_host2device();
    nonbond_a2.copy_host2device();

    // Setup cuda random
    // A state for each thread
    const int num_threads = batch_ligands * BLOCK_SIZE;
    curand_states.alloc(num_threads);

    // Simulate the population evolution for the given amount of time
    const auto num_generations = configuration.num_generations;
    // The shared memory contains:
    // - double scratch memory for intra block reductions
    // - each chromosome's score at last population evaluation

    call_kernel(batch_ligands,
                num_generations,
                configuration.tournament_length,
                configuration.mutation_prob,
                configuration.population_number,
                population_stride,
                batch_atoms,
                batch_rotamers,
                max_non_bonds,
                original_ligand_x.dev_pointer(),
                original_ligand_y.dev_pointer(),
                original_ligand_z.dev_pointer(),
                scratch_ligand_x.dev_pointer(),
                scratch_ligand_y.dev_pointer(),
                scratch_ligand_z.dev_pointer(),
                ligand_vol.dev_pointer(),
                ligand_solpar.dev_pointer(),
                ligand_charge.dev_pointer(),
                ligand_num_hbond.dev_pointer(),
                ligand_Rij_hb.dev_pointer(),
                ligand_Rii.dev_pointer(),
                ligand_epsij_hb.dev_pointer(),
                ligand_epsii.dev_pointer(),
                num_nonbonds.dev_pointer(),
                nonbond_a1.dev_pointer(),
                nonbond_a2.dev_pointer(),
                ligand_num_atoms.dev_pointer(),
                ligand_num_rotamers.dev_pointer(),
                ligand_fragments.dev_pointer(),
                frag_start_atom_indices.dev_pointer(),
                frag_stop_atom_indices.dev_pointer(),
                chromosomes.dev_pointer(),
                minimum_coord.x,
                minimum_coord.y,
                minimum_coord.z,
                maximum_coord.x,
                maximum_coord.y,
                maximum_coord.z,
                center.x,
                center.y,
                center.z,
                index_map.size_x(),
                index_map.size_xy(),
                atom_texs.wrappers_pointer.dev_pointer(),
                map_texture_index.dev_pointer(),
                electro_tex.dev_pointer(),
                desolv_tex.dev_pointer(),
                curand_states.dev_pointer(),
                ligand_scores.dev_pointer(),
                best_chromosomes.dev_pointer());

    // Copy back chromosomes and scores
    best_chromosomes.copy_device2host();
    ligand_scores.copy_device2host();

    // update the ligand position with the best one that we found
    index = 0;
    for (auto& ligand: std::span(incoming_batch.molecules.data(), incoming_batch.num_ligands)) {
      // Reset the random number generator to improve consistency
      apply(ligand.get()->get_x(),
            ligand.get()->get_y(),
            ligand.get()->get_z(),
            *(best_chromosomes.host_pointer() + index),
            batch_fragments[index]);
      ligand->properties.assign(property_type::SCORE, std::to_string(ligand_scores.host_pointer()[index]));
      ++index;
    }
  }
} // namespace mudock
