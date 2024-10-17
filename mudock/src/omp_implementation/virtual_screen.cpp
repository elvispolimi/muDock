#include <cstddef>
#include <cstring>
#include <cuda_runtime.h>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/grid.hpp>
#include <mudock/omp_implementation/evaluate_fitness.hpp>
#include <mudock/omp_implementation/virtual_screen.hpp>
#include <mudock/utils.hpp>
#include <omp.h>
#include <span>

namespace mudock {
  static constexpr std::size_t max_non_bonds{1024 * 10};

  void init_texture_memory(const grid_map &map, omp_object<fp_type> &tex_obj) {
    tex_obj.alloc(map.index.size_x() * map.index.size_y() * map.index.size_z());
    tex_obj.copy_host2device(map.data());
  }

  virtual_screen_omp::virtual_screen_omp(const knobs k,
                                         std::shared_ptr<const grid_atom_mapper> &grid_atom_maps,
                                         std::shared_ptr<const grid_map> &electro_map,
                                         std::shared_ptr<const grid_map> &desolv_map)
      : configuration(k),
        center(electro_map.get()->center),
        minimum_coord(electro_map.get()->minimum_coord),
        maximum_coord(electro_map.get()->maximum_coord),
        index_map(electro_map.get()->index) {
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
      omp_object<fp_type> &atom_tex = atom_texs.wrappers.back();
      const grid_map &grid_atom =
          grid_atom_maps.get()->get_atom_map(autodock_type_from_map(static_cast<device_map_textures>(index)));
      init_texture_memory(grid_atom, atom_tex);
      atom_texs.wrappers_pointer.host[index] = atom_tex.dev_pointer();
    }
    atom_texs.wrappers_pointer.copy_host2device();
  }

  void virtual_screen_omp::operator()(batch &incoming_batch) {
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
    const std::size_t total_chromosomes = configuration.population_number * batch_ligands;

    // Allocate memory for input data
    ligand_num_atoms.alloc(batch_ligands);
    ligand_num_rotamers.alloc(batch_ligands);
    original_ligand_x.alloc(tot_atoms_in_batch);
    original_ligand_y.alloc(tot_atoms_in_batch);
    original_ligand_z.alloc(tot_atoms_in_batch);
    ligand_vol.alloc(tot_atoms_in_batch);
    ligand_solpar.alloc(tot_atoms_in_batch);
    ligand_charge.alloc(tot_atoms_in_batch);
    ligand_Rij_hb.alloc(tot_atoms_in_batch);
    ligand_Rii.alloc(tot_atoms_in_batch);
    ligand_epsij_hb.alloc(tot_atoms_in_batch);
    ligand_epsii.alloc(tot_atoms_in_batch);
    ligand_num_hbond.alloc(tot_atoms_in_batch);
    // Bonds
    num_nonbonds.alloc(batch_ligands);
    nonbond_a1.alloc(batch_nonbonds);
    nonbond_a2.alloc(batch_nonbonds);
    // Fragments
    ligand_fragments.alloc(tot_rotamers_atoms_in_batch);
    frag_start_atom_indices.alloc(tot_rotamers_in_batch);
    frag_stop_atom_indices.alloc(tot_rotamers_in_batch);
    // Chromosomes and scores
    chromosomes.alloc(population_stride * batch_ligands);
    ligand_scores.alloc(batch_ligands);
    best_chromosomes.alloc(batch_ligands);
    // Scratchpad
    scratch_ligand_x.alloc(tot_atoms_in_batch);
    scratch_ligand_y.alloc(tot_atoms_in_batch);
    scratch_ligand_z.alloc(tot_atoms_in_batch);
    scratch_chromosome.alloc(total_chromosomes);
    // Support data precomputation
    map_texture_index.alloc(tot_atoms_in_batch);

    // Copy data
    std::size_t index{0};
    // TODO OpenMP NOTE not saved
    // Keep the fragments for the output
    // std::vector<fragments<static_containers>> batch_fragments;
    // batch_fragments.reserve(batch_ligands);
    // Support data structures
    grid<uint_fast8_t, index2D> nbmatrix{{static_cast<int>(batch_atoms), static_cast<int>(batch_atoms)}};
    std::vector<non_bond_parameter> non_bond_list;
    non_bond_list.reserve(max_non_bonds);
    for (auto &ligand: std::span(incoming_batch.molecules.data(), batch_ligands)) {
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

      std::memcpy((void *) (original_ligand_x.host_pointer() + stride_atoms),
                  x.data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (original_ligand_y.host_pointer() + stride_atoms),
                  y.data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (original_ligand_z.host_pointer() + stride_atoms),
                  z.data(),
                  num_atoms * sizeof(fp_type));
      // Fragments
      // Find out the rotatable bonds in the ligand
      auto graph = make_graph(ligand.get()->get_bonds());
      // TODO check this assignment
      const fragments<static_containers> l_fragments{graph,
                                                     ligand.get()->get_bonds(),
                                                     ligand.get()->num_atoms()};

      // Randomly initialize the population
      const auto num_rotamers                   = l_fragments.get_num_rotatable_bonds();
      ligand_num_rotamers.host_pointer()[index] = num_rotamers;
      const int stride_masks                    = index * batch_rotamers * batch_atoms;
      const int stride_rotamers                 = index * batch_rotamers;
      assert(static_cast<int>(batch_rotamers) > ligand.get()->num_rotamers());
      for (int rot = 0; rot < num_rotamers; ++rot) {
        std::memcpy((void *) (ligand_fragments.host_pointer() + stride_masks + rot * batch_atoms),
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
      if constexpr (is_debug())
        if (non_bond_list.size() >= max_non_bonds) {
          throw std::runtime_error("Bond list size exceed maximum value " +
                                   std::to_string(non_bond_list.size()) + ".");
        }

      num_nonbonds.host_pointer()[index] = non_bond_list.size();
      const int stride_nonbonds          = index * max_non_bonds;
      int nonbond_index{0};
      for (auto &bond: non_bond_list) {
        nonbond_a1.host_pointer()[stride_nonbonds + nonbond_index] = bond.a1;
        nonbond_a2.host_pointer()[stride_nonbonds + nonbond_index] = bond.a2;
        ++nonbond_index;
      }

      // Autodock typing
      std::memcpy((void *) (ligand_vol.host_pointer() + stride_atoms),
                  ligand.get()->get_vol().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_solpar.host_pointer() + stride_atoms),
                  ligand.get()->get_solpar().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_charge.host_pointer() + stride_atoms),
                  ligand.get()->get_charge().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_Rij_hb.host_pointer() + stride_atoms),
                  ligand.get()->get_Rij_hb().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_Rii.host_pointer() + stride_atoms),
                  ligand.get()->get_Rii().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_epsij_hb.host_pointer() + stride_atoms),
                  ligand.get()->get_epsij_hb().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_epsii.host_pointer() + stride_atoms),
                  ligand.get()->get_epsii().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_num_hbond.host_pointer() + stride_atoms),
                  ligand.get()->get_num_hbond().data(),
                  num_atoms * sizeof(int));

      std::size_t atom_index{0};
      for (auto &autodock_t: ligand.get()->get_autodock_type()) {
        const auto map_index = static_cast<int>(map_from_autodock_type(autodock_t));
        map_texture_index.host_pointer()[stride_atoms + atom_index] = map_index;
        ++atom_index;
      }

      ++index;
    }

    // Copy in
    ligand_num_atoms.copy_host2device();
    ligand_num_rotamers.copy_host2device();
    original_ligand_x.copy_host2device();
    original_ligand_y.copy_host2device();
    original_ligand_z.copy_host2device();
    ligand_vol.copy_host2device();
    ligand_solpar.copy_host2device();
    ligand_charge.copy_host2device();
    ligand_Rij_hb.copy_host2device();
    ligand_Rii.copy_host2device();
    ligand_epsij_hb.copy_host2device();
    ligand_epsii.copy_host2device();
    ligand_num_hbond.copy_host2device();
    num_nonbonds.copy_host2device();
    nonbond_a1.copy_host2device();
    nonbond_a2.copy_host2device();
    ligand_fragments.copy_host2device();
    frag_start_atom_indices.copy_host2device();
    frag_stop_atom_indices.copy_host2device();
    map_texture_index.copy_host2device();

    // Setup OpenMP random data
    // A state for each chromosome
    omp_states.alloc(configuration.population_number * batch_ligands);

    // Get device pointers
    const auto *d_num_atoms               = ligand_num_atoms.dev_pointer();
    const auto *d_num_rotamers            = ligand_num_rotamers.dev_pointer();
    const auto *d_original_ligand_x       = original_ligand_x.dev_pointer();
    const auto *d_original_ligand_y       = original_ligand_y.dev_pointer();
    const auto *d_original_ligand_z       = original_ligand_z.dev_pointer();
    auto *d_scratch_ligand_x              = scratch_ligand_x.dev_pointer();
    auto *d_scratch_ligand_y              = scratch_ligand_y.dev_pointer();
    auto *d_scratch_ligand_z              = scratch_ligand_z.dev_pointer();
    auto *d_scratch_chromosome            = scratch_chromosome.dev_pointer();
    const auto *d_ligand_vol              = ligand_vol.dev_pointer();
    const auto *d_ligand_solpar           = ligand_solpar.dev_pointer();
    const auto *d_ligand_charge           = ligand_charge.dev_pointer();
    const auto *d_ligand_Rij_hb           = ligand_Rij_hb.dev_pointer();
    const auto *d_ligand_Rii              = ligand_Rii.dev_pointer();
    const auto *d_ligand_epsij_hb         = ligand_epsij_hb.dev_pointer();
    const auto *d_ligand_epsii            = ligand_epsii.dev_pointer();
    const auto *d_ligand_num_hbond        = ligand_num_hbond.dev_pointer();
    const auto *d_num_nonbonds            = num_nonbonds.dev_pointer();
    const auto *d_nonbond_a1              = nonbond_a1.dev_pointer();
    const auto *d_nonbond_a2              = nonbond_a2.dev_pointer();
    const auto *d_ligand_fragments        = ligand_fragments.dev_pointer();
    const auto *d_frag_start_atom_indices = frag_start_atom_indices.dev_pointer();
    const auto *d_frag_stop_atom_indices  = frag_stop_atom_indices.dev_pointer();
    auto *d_chromosomes                   = chromosomes.dev_pointer();
    const auto *d_atom_textures           = atom_texs.wrappers_pointer.dev_pointer();
    const auto *d_atom_tex_indexes        = map_texture_index.dev_pointer();
    const auto *d_electro_texture         = electro_tex.dev_pointer();
    const auto *d_desolv_texture          = desolv_tex.dev_pointer();
    auto *d_omp_states                    = omp_states.dev_pointer();
    auto *d_ligand_scores                 = ligand_scores.dev_pointer();
    auto *d_best_chromosomes              = best_chromosomes.dev_pointer();

    // Execute the kernel
    evaluate_fitness(batch_ligands,
                     configuration.num_generations,
                     configuration.tournament_length,
                     configuration.mutation_prob,
                     configuration.population_number,
                     population_stride,
                     batch_atoms,
                     batch_rotamers,
                     max_non_bonds,
                     d_original_ligand_x,
                     d_original_ligand_y,
                     d_original_ligand_z,
                     d_scratch_ligand_x,
                     d_scratch_ligand_y,
                     d_scratch_ligand_z,
                     d_scratch_chromosome,
                     d_ligand_vol,
                     d_ligand_solpar,
                     d_ligand_charge,
                     d_ligand_num_hbond,
                     d_ligand_Rij_hb,
                     d_ligand_Rii,
                     d_ligand_epsij_hb,
                     d_ligand_epsii,
                     d_num_nonbonds,
                     d_nonbond_a1,
                     d_nonbond_a2,
                     d_num_atoms,
                     d_num_rotamers,
                     d_ligand_fragments,
                     d_frag_start_atom_indices,
                     d_frag_stop_atom_indices,
                     d_chromosomes,
                     minimum_coord,
                     maximum_coord,
                     center,
                     index_map,
                     d_atom_textures,
                     d_atom_tex_indexes,
                     d_electro_texture,
                     d_desolv_texture,
                     d_omp_states,
                     d_ligand_scores,
                     d_best_chromosomes);

    // Copy back chromosomes and scores
    best_chromosomes.copy_device2host();
    ligand_scores.copy_device2host();

    // update the ligand position with the best one that we found
    index = 0;
    for (auto &ligand: std::span(incoming_batch.molecules.data(), incoming_batch.num_ligands)) {
      auto graph = make_graph(ligand.get()->get_bonds());
      // TODO check this assignment
      const fragments<static_containers> l_fragments{graph,
                                                     ligand.get()->get_bonds(),
                                                     ligand.get()->num_atoms()};
      // Reset the random number generator to improve consistency
      apply(ligand.get()->get_x(),
            ligand.get()->get_y(),
            ligand.get()->get_z(),
            *(best_chromosomes.host_pointer() + index),
            l_fragments);
      ligand->properties.assign(property_type::SCORE, std::to_string(ligand_scores.host_pointer()[index]));
      ++index;
    }
  } // namespace mudock
} // namespace mudock
