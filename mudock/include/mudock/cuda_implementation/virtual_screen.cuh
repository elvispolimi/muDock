#pragma once

#include <mudock/batch.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cuda_implementation/cuda_random.cuh>
#include <mudock/cuda_implementation/cuda_wrapper.cuh>
#include <mudock/cuda_implementation/device.cuh>
#include <mudock/grid.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <vector>

namespace mudock {

  class virtual_screen_cuda {
    // the configuration of the GA algorithm
    knobs configuration;

    std::shared_ptr<const device> dev;
    cudaStream_wrapper stream;

    // Data area
    // TODO some of these can be placed into shared memory
    cuda_wrapper<std::vector, fp_type> original_ligand_x, original_ligand_y, original_ligand_z,
        scratch_ligand_x, scratch_ligand_y, scratch_ligand_z, ligand_vol, ligand_solpar, ligand_charge;
    cuda_wrapper<std::vector, int> ligand_num_atoms, ligand_num_rotamers;
    // Fragments
    cuda_wrapper<std::vector, int> ligand_fragments, ligand_fragments_start;
    cuda_wrapper<std::vector, int> frag_start_atom_indices, frag_stop_atom_indices, frag_indices_start;
    // Non-bonds
    cuda_wrapper<std::vector, int> index_nonbonds, nonbond_a1, nonbond_a2, nonbond_xB;
    cuda_wrapper<std::vector, fp_type> nonbond_cA, nonbond_cB;

    // CUDA data precomputation
    cuda_wrapper<std::vector, int> map_texture_index;

    // Return energy
    cuda_wrapper<std::vector, fp_type> ligand_scores;

    // define the GA population
    cuda_wrapper<std::vector, chromosome> chromosomes;
    cuda_wrapper<std::vector, chromosome> best_chromosomes;

    // Random generation
    cuda_random_object curand_states;

  public:
    virtual_screen_cuda(const knobs k, const std::shared_ptr<const device> dev);

    void operator()(batch& incoming_batch);
  };
} // namespace mudock
