#pragma once

#include <mudock/batch.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cuda_implementation/map_textures.cuh>
#include <mudock/grid.hpp>
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

    // Data area
    // TODO some of these can be placed into shared memory
    hip_wrapper<std::vector, fp_type> original_ligand_x, original_ligand_y, original_ligand_z,
        scratch_ligand_x, scratch_ligand_y, scratch_ligand_z, ligand_vol, ligand_solpar, ligand_charge,
        ligand_Rij_hb, ligand_Rii, ligand_epsij_hb, ligand_epsii;
    hip_wrapper<std::vector, int> ligand_num_hbond, ligand_num_atoms, ligand_num_rotamers;
    // Fragments
    hip_wrapper<std::vector, int> ligand_fragments;
    hip_wrapper<std::vector, int> frag_start_atom_indices, frag_stop_atom_indices;
    // Non-onds
    hip_wrapper<std::vector, int> num_nonbonds, nonbond_a1, nonbond_a2;

    // HIP data precomputation
    hip_wrapper<std::vector, int> map_texture_index;

    // Return energy
    hip_wrapper<std::vector, fp_type> ligand_scores;

    // define the GA population
    hip_wrapper<std::vector, chromosome> chromosomes;
    hip_wrapper<std::vector, chromosome> best_chromosomes;

    // Grid Maps
    const point3D center_maps;
    hipTextureObject_t electro_tex, desolv_tex;
    hip_wrapper<std::vector, hipTextureObject_t> atom_texs;

    // Random generation
    hip_random_object curand_states;

    // Get the wavefront size of the first device
    int wavefront_size;

  public:
    virtual_screen_hip(const knobs k,
                       std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                       std::shared_ptr<const grid_map>& electro_map,
                       std::shared_ptr<const grid_map>& desolv_map,
                       const int device_id);

    void operator()(batch& incoming_batch);
  };
} // namespace mudock
