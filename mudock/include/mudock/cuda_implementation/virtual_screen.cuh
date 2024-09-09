#pragma once

#include <cstddef>
#include <memory>
#include <mudock/batch.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cuda_implementation/cuda_object.cuh>
#include <mudock/cuda_implementation/cuda_wrapper.cuh>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <random>
#include <vector>

namespace mudock {

  class virtual_screen_cuda {
    // this will be the scratchpad memory for the CUDA implementation
    // NOTE: we will implement a random scoring function to test the infrastructure
    std::mt19937 generator;
    std::uniform_real_distribution<fp_type> dist;

    // the configuration of the GA algorithm
    knobs configuration;

    // Data area
    // TODO some of these can be placed into shared memory
    cuda_wrapper<std::vector, fp_type> ligand_x, ligand_y, ligand_z, ligand_vol, ligand_solpar, ligand_charge,
        ligand_Rij_hb, ligand_Rii, ligand_epsij_hb, ligand_epsii;
    cuda_wrapper<std::vector, std::size_t> ligand_autodock_type, ligand_num_hbond, ligand_num_atoms, ligand_num_rotamers;
    cuda_wrapper<std::vector, typename fragments<static_containers>::value_type> ligand_fragments;
    cuda_wrapper<std::vector, std::size_t> frag_start_atom_indices, frag_stop_atom_indices;
    // TODO maps

    // Return energy
    cuda_wrapper<std::vector, fp_type> ligand_scores;

    // define the GA population
    cuda_wrapper<std::vector, individual> population, next_population;

  public:
    virtual_screen_cuda(const knobs k);

    void operator()(batch& incoming_batch);
  };
} // namespace mudock
