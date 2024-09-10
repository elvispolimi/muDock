#include <alloca.h>
#include <cstddef>
#include <cstring>
#include <cuda_runtime.h>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cuda_implementation/evaluate_fitness.cuh>
#include <mudock/cuda_implementation/virtual_screen.cuh>
#include <mudock/grid.hpp>
#include <span>

#define BLOCK_SIZE 32

namespace mudock {

  void init_texture_memory(std::shared_ptr<const grid_map> &map, cudaTextureObject_t tex_obj) {
    // Create 3D CUDA array for the texture
    cudaArray *d_array;
    cudaChannelFormatDesc channel_desc = cudaCreateChannelDesc<fp_type>();

    const index3D map_index = map.get()->index;
    cudaExtent extent       = make_cudaExtent(map_index.size_x(), map_index.size_y(), map_index.size_z());
    cudaMalloc3DArray(&d_array, &channel_desc, extent);

    // Copy data from host to the 3D CUDA array
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr            = make_cudaPitchedPtr((void *) map.get()->data(),
                                            map_index.size_x() * sizeof(fp_type),
                                            map_index.size_y(),
                                            map_index.size_z());
    copyParams.dstArray          = d_array;
    copyParams.extent            = extent;
    copyParams.kind              = cudaMemcpyHostToDevice;
    cudaMemcpy3D(&copyParams);

    // Create texture object
    struct cudaResourceDesc res_desc;
    memset(&res_desc, 0, sizeof(res_desc));
    res_desc.resType         = cudaResourceTypeArray;
    res_desc.res.array.array = d_array;

    struct cudaTextureDesc tex_desc;
    memset(&tex_desc, 0, sizeof(tex_desc));
    tex_desc.addressMode[0]   = cudaAddressModeClamp;
    tex_desc.addressMode[1]   = cudaAddressModeClamp;
    tex_desc.addressMode[2]   = cudaAddressModeClamp;
    tex_desc.filterMode       = cudaFilterModeLinear; // Enable linear interpolation
    tex_desc.readMode         = cudaReadModeElementType;
    tex_desc.normalizedCoords = 0; // We will use unnormalized coordinates

    // Create the texture object
    cudaCreateTextureObject(&tex_obj, &res_desc, &tex_desc, NULL);
  }

  virtual_screen_cuda::virtual_screen_cuda(const knobs k,
                                           std::shared_ptr<const grid_atom_mapper> &grid_atom_maps,
                                           std::shared_ptr<const grid_map> &electro_map,
                                           std::shared_ptr<const grid_map> &desolv_map)
      : dist(fp_type{0}, fp_type{10}), configuration(k) {
    // Allocate grid maps
    init_texture_memory(electro_map, electro_tex);
  }

  void virtual_screen_cuda::operator()(batch &incoming_batch) {
    const std::size_t batch_atoms    = incoming_batch.batch_max_atoms;
    const std::size_t batch_rotamers = incoming_batch.batch_max_rotamers;
    const std::size_t batch_ligands  = incoming_batch.num_ligands;
    // Resize data structures
    const std::size_t tot_atoms_in_batch          = batch_ligands * batch_atoms;
    const std::size_t tot_atoms_in_population     = tot_atoms_in_batch * configuration.population_number;
    const std::size_t tot_rotamers_atoms_in_batch = tot_atoms_in_batch * batch_rotamers;
    const std::size_t tot_rotamers_in_batch       = batch_ligands * batch_rotamers;
    ligand_x.alloc(tot_atoms_in_population, fp_type{0});
    ligand_y.alloc(tot_atoms_in_population, fp_type{0});
    ligand_z.alloc(tot_atoms_in_population, fp_type{0});
    ligand_fragments.alloc(tot_rotamers_atoms_in_batch);
    frag_start_atom_indices.alloc(tot_rotamers_in_batch);
    frag_stop_atom_indices.alloc(tot_rotamers_in_batch);
    ligand_vol.alloc(tot_atoms_in_batch, fp_type{0});
    ligand_solpar.alloc(tot_atoms_in_batch, fp_type{0});
    ligand_charge.alloc(tot_atoms_in_batch, fp_type{0});
    ligand_Rij_hb.alloc(tot_atoms_in_batch, fp_type{0});
    ligand_Rii.alloc(tot_atoms_in_batch, fp_type{0});
    ligand_epsij_hb.alloc(tot_atoms_in_batch, fp_type{0});
    ligand_epsii.alloc(tot_atoms_in_batch, fp_type{0});
    ligand_autodock_type.alloc(tot_atoms_in_batch, static_cast<std::size_t>(autodock_ff::H));
    ligand_num_hbond.alloc(batch_ligands, std::size_t{0});
    ligand_num_atoms.alloc(batch_ligands, std::size_t{0});
    ligand_num_rotamers.alloc(batch_ligands, std::size_t{0});
    ligand_scores.alloc(batch_ligands, fp_type{0});
    // GA Data structures
    population.alloc(configuration.population_number * batch_ligands);

    // Define the range that we can use to mutate and generate the chromosome
    auto init_change_distribution     = std::uniform_int_distribution(-45, 45);
    auto mutation_change_distribution = std::uniform_int_distribution(-10, 10);
    auto mutation_coin_distribution   = std::uniform_real_distribution{fp_type{0}, fp_type{1.0}};
    const auto coordinate_step        = fp_type{0.2};
    const auto angle_step             = fp_type{4};

    // Copy data
    std::size_t index{0};
    for (auto &ligand: std::span(incoming_batch.molecules.data(), batch_ligands)) {
      // Atoms and bonds
      const std::size_t num_atoms = ligand.get()->num_atoms();
      std::memcpy((void *) (ligand_num_atoms.host_pointer() + index), &num_atoms, sizeof(std::size_t));
      // TODO bonds
      // Place the molecule to the center of the target protein
      const auto x = ligand.get()->get_x(), y = ligand.get()->get_y(), z = ligand.get()->get_z();
      const auto ligand_center_of_mass = compute_center_of_mass(x, y, z);
      // TODO missing
      // translate_molecule(x,
      //                    y,
      //                    z,
      //                    electro_map->center.x - ligand_center_of_mass.x,
      //                    electro_map->center.y - ligand_center_of_mass.y,
      //                    electro_map->center.z - ligand_center_of_mass.z);
      // Coordinates
      for (std::size_t i = 0; i < configuration.population_number; ++i) {
        std::memcpy((void *) (ligand_x.host_pointer() +
                              index * batch_atoms * configuration.population_number + i * batch_atoms),
                    ligand.get()->get_x().data(),
                    num_atoms * sizeof(fp_type));
        std::memcpy((void *) (ligand_y.host_pointer() +
                              index * batch_atoms * configuration.population_number + i * batch_atoms),
                    ligand.get()->get_y().data(),
                    num_atoms * sizeof(fp_type));
        std::memcpy((void *) (ligand_z.host_pointer() +
                              index * batch_atoms * configuration.population_number + i * batch_atoms),
                    ligand.get()->get_z().data(),
                    num_atoms * sizeof(fp_type));
      }
      // Fragments
      // Find out the rotatable bonds in the ligand
      auto graph = make_graph(ligand.get()->get_bonds());
      const fragments<static_containers> l_fragments{graph,
                                                     ligand.get()->get_bonds(),
                                                     ligand.get()->num_atoms()};

      // Randomly initialize the population
      const auto num_rotamers = l_fragments.get_num_rotatable_bonds();
      std::memcpy((void *) (ligand_num_rotamers.host_pointer() + index), &num_rotamers, sizeof(std::size_t));
      for (std::size_t rot = 0; rot < num_rotamers; ++rot) {
        std::memcpy((void *) (ligand_fragments.host_pointer() + index * batch_rotamers * batch_atoms +
                              rot * batch_atoms),
                    l_fragments.get_mask(rot).data(),
                    num_atoms * sizeof(fp_type));
        const auto [start_index, stop_index] = l_fragments.get_rotatable_atoms(rot);
        std::memcpy((void *) (frag_start_atom_indices.host_pointer() + index * batch_rotamers),
                    &start_index,
                    sizeof(std::size_t));
        std::memcpy((void *) (frag_stop_atom_indices.host_pointer() + index * batch_rotamers),
                    &stop_index,
                    sizeof(std::size_t));
      }
      // Autodock typing
      std::memcpy((void *) (ligand_vol.host_pointer() + index * batch_atoms),
                  ligand.get()->get_vol().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_solpar.host_pointer() + index * batch_atoms),
                  ligand.get()->get_solpar().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_charge.host_pointer() + index * batch_atoms),
                  ligand.get()->get_charge().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_Rij_hb.host_pointer() + index * batch_atoms),
                  ligand.get()->get_Rij_hb().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_Rii.host_pointer() + index * batch_atoms),
                  ligand.get()->get_Rii().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_epsij_hb.host_pointer() + index * batch_atoms),
                  ligand.get()->get_epsij_hb().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_epsii.host_pointer() + index * batch_atoms),
                  ligand.get()->get_epsii().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_autodock_type.host_pointer() + index * batch_atoms),
                  ligand.get()->get_autodock_type().data(),
                  num_atoms * sizeof(std::size_t));
      std::memcpy((void *) (ligand_num_hbond.host_pointer() + index),
                  ligand.get()->get_num_hbond().data(),
                  num_atoms * sizeof(std::size_t));

      for (std::size_t i = 0; i < configuration.population_number; ++i) {
        auto &element = population.host[index * configuration.population_number + i];
        for (std::size_t i{0}; i < 3; ++i) { // initialize the rigid translation
          element.genes[i] = static_cast<fp_type>(init_change_distribution(generator)) * coordinate_step;
        }
        for (std::size_t i{3}; i < 3 + num_rotamers; ++i) { // initialize the rotations
          element.genes[i] = static_cast<fp_type>(init_change_distribution(generator)) * angle_step;
        }
      }

      ++index;
    }
    // Copy in
    ligand_num_atoms.copy_host2device();
    ligand_x.copy_host2device();
    ligand_y.copy_host2device();
    ligand_z.copy_host2device();
    ligand_fragments.copy_host2device();
    frag_start_atom_indices.copy_host2device();
    frag_stop_atom_indices.copy_host2device();
    ligand_vol.copy_host2device();
    ligand_solpar.copy_host2device();
    ligand_charge.copy_host2device();
    ligand_Rij_hb.copy_host2device();
    ligand_Rii.copy_host2device();
    ligand_epsij_hb.copy_host2device();
    ligand_epsii.copy_host2device();
    ligand_autodock_type.copy_host2device();
    ligand_num_hbond.copy_host2device();
    population.copy_host2device();

    // Simulate the population evolution for the given amount of time
    const auto num_generations = configuration.num_generations;
    for (std::size_t generation = 0; generation < num_generations; ++generation) {
      // TODO checks if everything fit into shared memory
      evaluate_fitness<<<batch_ligands, BLOCK_SIZE>>>(configuration.population_number,
                                                      batch_atoms,
                                                      batch_rotamers,
                                                      ligand_x.dev_pointer(),
                                                      ligand_y.dev_pointer(),
                                                      ligand_z.dev_pointer(),
                                                      ligand_vol.dev_pointer(),
                                                      ligand_solpar.dev_pointer(),
                                                      ligand_charge.dev_pointer(),
                                                      ligand_num_hbond.dev_pointer(),
                                                      ligand_Rij_hb.dev_pointer(),
                                                      ligand_Rii.dev_pointer(),
                                                      ligand_epsij_hb.dev_pointer(),
                                                      ligand_epsii.dev_pointer(),
                                                      ligand_autodock_type.dev_pointer(),
                                                      NULL,
                                                      ligand_num_atoms.dev_pointer(),
                                                      ligand_num_rotamers.dev_pointer(),
                                                      ligand_fragments.dev_pointer(),
                                                      frag_start_atom_indices.dev_pointer(),
                                                      frag_stop_atom_indices.dev_pointer(),
                                                      population.dev_pointer(),
                                                      NULL,
                                                      NULL,
                                                      NULL,
                                                      ligand_scores.dev_pointer());
      MUDOCK_CHECK_KERNELCALL();
      MUDOCK_CHECK(cudaDeviceSynchronize());

      // Generate the new population
      // TODO
    }

    // update the ligand position with the best one that we found
    // TODO

    // for (auto& ligand: std::span(incoming_batch.molecules.data(), incoming_batch.num_ligands)) {
    //   // Reset the random number generator to improve consistency
    //   generator = std::mt19937{ligand.get()->num_atoms()};
    //   ligand->properties.assign(property_type::SCORE, std::to_string(dist(generator)));
    // }
  }
} // namespace mudock
