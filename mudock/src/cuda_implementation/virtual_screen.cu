#include <alloca.h>
#include <cstddef>
#include <cstring>
#include <cuda_runtime.h>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/cuda_implementation/evaluate_fitness.cuh>
#include <mudock/cuda_implementation/virtual_screen.cuh>
#include <mudock/grid.hpp>
#include <mudock/utils.hpp>
#include <span>

// Keep it to 32 to enable warp optimizations
#define BLOCK_SIZE 32
namespace mudock {
  static constexpr std::size_t max_non_bonds{1024 * 10};

  void init_texture_memory(const grid_map &map, cudaTextureObject_t &tex_obj) {
    // Create 3D CUDA array for the texture
    cudaArray *d_array;
    cudaChannelFormatDesc channel_desc = cudaCreateChannelDesc<fp_type>();

    const index3D map_index = map.index;
    cudaExtent extent       = make_cudaExtent(map_index.size_x(), map_index.size_y(), map_index.size_z());
    MUDOCK_CHECK(cudaMalloc3DArray(&d_array, &channel_desc, extent));

    // Copy data from host to the 3D CUDA array
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr =
        make_cudaPitchedPtr((void *) map.data(), extent.width * sizeof(fp_type), extent.width, extent.height);
    copyParams.srcArray = nullptr;
    copyParams.srcPos   = make_cudaPos(0, 0, 0);
    copyParams.dstArray = d_array;
    copyParams.dstPtr   = cudaPitchedPtr{};
    copyParams.dstPos   = make_cudaPos(0, 0, 0);
    copyParams.extent   = extent;
    copyParams.kind     = cudaMemcpyHostToDevice;
    MUDOCK_CHECK(cudaMemcpy3D(&copyParams));

    // Create texture object
    cudaResourceDesc res_desc;
    memset(&res_desc, 0, sizeof(cudaResourceDesc));
    res_desc.resType         = cudaResourceTypeArray;
    res_desc.res.array.array = d_array;

    cudaTextureDesc tex_desc;
    memset(&tex_desc, 0, sizeof(cudaTextureDesc));
    tex_desc.addressMode[0] = cudaAddressModeClamp;
    tex_desc.addressMode[1] = cudaAddressModeClamp;
    tex_desc.addressMode[2] = cudaAddressModeClamp;
    // tex_desc.filterMode     = cudaFilterModeLinear; // Enable linear interpolation
    tex_desc.filterMode       = cudaFilterModePoint;
    tex_desc.readMode         = cudaReadModeElementType;
    tex_desc.normalizedCoords = false; // We will use unnormalized coordinates

    // Create the texture object
    MUDOCK_CHECK(cudaCreateTextureObject(&tex_obj, &res_desc, &tex_desc, NULL));
  }

  virtual_screen_cuda::virtual_screen_cuda(const knobs k,
                                           std::shared_ptr<const grid_atom_mapper> &grid_atom_maps,
                                           std::shared_ptr<const grid_map> &electro_map,
                                           std::shared_ptr<const grid_map> &desolv_map)
      : configuration(k), center_maps(electro_map.get()->center) {
    // TODO move this part into the cuda_worker -> once per GPU
    // Allocate grid maps
    assert(electro_map.get()->index == desolv_map.get()->index &&
           desolv_map.get()->index == grid_atom_maps.get()->get_index());
    init_texture_memory(*electro_map.get(), electro_tex);

    init_texture_memory(*desolv_map.get(), desolv_tex);

    std::size_t index{0};
    atom_texs.alloc(num_device_map_textures());
    for (auto &atom_tex: atom_texs.host) {
      const grid_map &grid_atom =
          grid_atom_maps.get()->get_atom_map(autodock_type_from_map(static_cast<device_map_textures>(index)));
      init_texture_memory(grid_atom, atom_tex);
      ++index;
    }
    atom_texs.copy_host2device();

    // Grid spacing fixed to 0.5 Angstrom
    setup_constant_memory(electro_map.get()->minimum_coord,
                          electro_map.get()->maximum_coord,
                          electro_map.get()->center);
  }

  void virtual_screen_cuda::operator()(batch &incoming_batch) {
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
                         center_maps.x - ligand_center_of_mass.x,
                         center_maps.y - ligand_center_of_mass.y,
                         center_maps.z - ligand_center_of_mass.z);

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
      batch_fragments[index]  = {graph, ligand.get()->get_bonds(), ligand.get()->num_atoms()};
      const auto &l_fragments = batch_fragments[index];

      // Randomly initialize the population
      const auto num_rotamers                   = l_fragments.get_num_rotatable_bonds();
      ligand_num_rotamers.host_pointer()[index] = num_rotamers;
      const int stride_masks                    = index * batch_rotamers * batch_atoms;
      const int stride_rotamers                 = index * batch_rotamers;
      assert(batch_rotamers > ligand.get()->num_rotamers());
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
    // TODO checks if everything fit into shared memory
    // The shared memory contains:
    // - each chromosome's score at last population evaluation
    // Max due to the reduction at the end to find the highest scores per each chromosome
    const std::size_t min_energy_reduction_s_mem =
        std::max(configuration.population_number, static_cast<std::size_t>(BLOCK_SIZE)) * sizeof(fp_type);
    const std::size_t shared_mem = min_energy_reduction_s_mem;
    evaluate_fitness<<<batch_ligands, BLOCK_SIZE, shared_mem>>>(num_generations,
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
                                                                atom_texs.dev_pointer(),
                                                                map_texture_index.dev_pointer(),
                                                                electro_tex,
                                                                desolv_tex,
                                                                curand_states.dev_pointer(),
                                                                ligand_scores.dev_pointer(),
                                                                best_chromosomes.dev_pointer());
    MUDOCK_CHECK_KERNELCALL();
    MUDOCK_CHECK(cudaDeviceSynchronize());

    // Copy back chromosomes and scores
    best_chromosomes.copy_device2host();
    ligand_scores.copy_device2host();

    // update the ligand position with the best one that we found
    index = 0;
    for (auto &ligand: std::span(incoming_batch.molecules.data(), incoming_batch.num_ligands)) {
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
