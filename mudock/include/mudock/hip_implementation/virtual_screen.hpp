#pragma once

#include <mudock/batch.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/grid.hpp>
#include <mudock/hip_implementation/device.hpp>
#include <mudock/hip_implementation/hip_random.hpp>
#include <mudock/hip_implementation/hip_wrapper.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <vector>

namespace mudock {

  class virtual_screen_hip {
    // the configuration of the GA algorithm
    knobs configuration;
    std::shared_ptr<const device> dev;
    hipStream_wrapper stream;

    // Data area
    // TODO some of these can be placed into shared memory
    hip_wrapper<std::vector, fp_type> original_ligand_x, original_ligand_y, original_ligand_z,
        scratch_ligand_x, scratch_ligand_y, scratch_ligand_z, ligand_vol, ligand_solpar, ligand_charge;
    hip_wrapper<std::vector, int> ligand_num_hbond, ligand_num_atoms, ligand_num_rotamers;
    // Fragments
    hip_wrapper<std::vector, int> ligand_fragments, ligand_fragments_start;
    hip_wrapper<std::vector, int> frag_start_atom_indices, frag_stop_atom_indices, frag_indices_start;
    // Non-bonds
    hip_wrapper<std::vector, int> index_nonbonds, nonbond_a1, nonbond_a2, nonbond_xB;
    hip_wrapper<std::vector, fp_type> nonbond_cA, nonbond_cB;

    // HIP data precomputation
    hip_wrapper<std::vector, int> map_texture_index;

    // Return energy
    hip_wrapper<std::vector, fp_type> ligand_scores;

    // define the GA population
    hip_wrapper<std::vector, chromosome> chromosomes;
    hip_wrapper<std::vector, chromosome> best_chromosomes;

    // Random generation
    // TODO check random generation on HIP performance
    hip_random_object hiprand_states;

  public:
    virtual_screen_hip(const knobs k, const std::shared_ptr<const device> dev);

    void operator()(batch<autodock_ligand>& incoming_batch);
  };
} // namespace mudock
