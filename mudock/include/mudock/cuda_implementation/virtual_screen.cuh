#pragma once

#include <cstddef>
#include <memory>
#include <mudock/batch.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cuda_implementation/cuda_object.cuh>
#include <mudock/cuda_implementation/cuda_wrapper.cuh>
#include <mudock/cuda_implementation/cuda_random.cuh>
#include <mudock/cuda_implementation/map_textures.cuh>
#include <mudock/grid.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <random>
#include <vector>

namespace mudock {

  class virtual_screen_cuda {
    // the configuration of the GA algorithm
    knobs configuration;

    // Data area
    // TODO some of these can be placed into shared memory
    cuda_wrapper<std::vector, fp_type> ligand_x, ligand_y, ligand_z, ligand_vol, ligand_solpar, ligand_charge,
        ligand_Rij_hb, ligand_Rii, ligand_epsij_hb, ligand_epsii;
    cuda_wrapper<std::vector, int> ligand_num_hbond, ligand_num_atoms, ligand_num_rotamers;
    // Fragments
    cuda_wrapper<std::vector, int> ligand_fragments;
    cuda_wrapper<std::vector, int> frag_start_atom_indices, frag_stop_atom_indices;
    // Non-onds
    cuda_wrapper<std::vector, int> num_nonbonds, nonbond_a1, nonbond_a2;

    // CUDA data precomputation
    cuda_wrapper<std::vector, int> map_texture_index;

    // Return energy
    cuda_wrapper<std::vector, fp_type> ligand_scores;

    // define the GA population
    cuda_wrapper<std::vector, chromosome> chromosomes;

    // Grid Maps
    const index3D index_maps;
    const point3D center_maps;
    cudaTextureObject_t electro_tex, desolv_tex;
    cuda_wrapper<std::vector, cudaTextureObject_t> atom_texs;

    // Random generation
    cuda_random_object curand_states;

  public:
    virtual_screen_cuda(const knobs k,
                        std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                        std::shared_ptr<const grid_map>& electro_map,
                        std::shared_ptr<const grid_map>& desolv_map);

    void operator()(batch& incoming_batch);
  };
} // namespace mudock
