#pragma once

#include <mudock/batch.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cuda_implementation/map_textures.cuh>
#include <mudock/grid.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <mudock/omp_implementation/omp_random.hpp>
#include <mudock/omp_implementation/omp_wrapper.hpp>
#include <mudock/type_alias.hpp>
#include <vector>

namespace mudock {
  struct wrappers_container {
    std::vector<omp_object<fp_type>> wrappers;
    omp_wrapper<std::vector, fp_type*> wrappers_pointer;

    wrappers_container(): wrappers_pointer(){};
  };

  class virtual_screen_omp {
    // the configuration of the GA algorithm
    knobs configuration;

    // Data area
    // TODO some of these can be placed into shared memory
    omp_wrapper<std::vector, fp_type> original_ligand_x, original_ligand_y, original_ligand_z, ligand_vol,
        ligand_solpar, ligand_charge, ligand_Rij_hb, ligand_Rii, ligand_epsij_hb, ligand_epsii;
    omp_wrapper<std::vector, int> ligand_num_hbond, ligand_num_atoms, ligand_num_rotamers;
    // Fragments
    omp_wrapper<std::vector, int> ligand_fragments;
    omp_wrapper<std::vector, int> frag_start_atom_indices, frag_stop_atom_indices;
    // Non-bonds
    omp_wrapper<std::vector, int> num_nonbonds, nonbond_a1, nonbond_a2;

    // Return energy
    omp_wrapper<std::vector, fp_type> ligand_scores;

    // Define the GA population
    omp_wrapper<std::vector, chromosome> chromosomes;
    omp_wrapper<std::vector, chromosome> best_chromosomes;

    // Grid Maps
    const point3D center, minimum_coord, maximum_coord;
    const index3D index_map;
    omp_object<fp_type> electro_tex, desolv_tex;
    wrappers_container atom_texs;
    omp_wrapper<std::vector, int> map_texture_index;

    // Scratchpads
    omp_object<fp_type> scratch_chromosome;
    omp_object<fp_type> scratch_ligand_x, scratch_ligand_y, scratch_ligand_z;

    // Random generation
    omp_random_object omp_states;

  public:
    virtual_screen_omp(const knobs k,
                       std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                       std::shared_ptr<const grid_map>& electro_map,
                       std::shared_ptr<const grid_map>& desolv_map);

    void operator()(batch& incoming_batch);
  };
} // namespace mudock
