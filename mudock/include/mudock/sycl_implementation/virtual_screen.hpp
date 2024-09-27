#pragma once

#include <mudock/batch.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/grid.hpp>
#include <mudock/knobs.hpp>
#include <mudock/sycl_implementation/sycl_random.hpp>
#include <mudock/sycl_implementation/sycl_wrapper.hpp>

namespace mudock {
  struct wrappers_container {
    std::vector<sycl_object<fp_type>> wrappers;
    sycl_wrapper<std::vector, fp_type*> wrappers_pointer;

    wrappers_container(sycl::queue& queue): wrappers_pointer(queue){};
  };

  class virtual_screen_sycl {
    // SYCL queue
    sycl::queue& queue;
    int subgroup_size;

    // the configuration of the GA algorithm
    knobs configuration;

    // Data area
    // TODO some of these can be placed into shared memory
    sycl_wrapper<std::vector, fp_type> original_ligand_x, original_ligand_y, original_ligand_z,
        scratch_ligand_x, scratch_ligand_y, scratch_ligand_z, ligand_vol, ligand_solpar, ligand_charge,
        ligand_Rij_hb, ligand_Rii, ligand_epsij_hb, ligand_epsii;
    sycl_wrapper<std::vector, int> ligand_num_hbond, ligand_num_atoms, ligand_num_rotamers;
    // Fragments
    sycl_wrapper<std::vector, int> ligand_fragments;
    sycl_wrapper<std::vector, int> frag_start_atom_indices, frag_stop_atom_indices;
    // Non-onds
    sycl_wrapper<std::vector, int> num_nonbonds, nonbond_a1, nonbond_a2;

    // SYCL data precomputation
    sycl_wrapper<std::vector, int> map_texture_index;

    // Return energy
    sycl_wrapper<std::vector, fp_type> ligand_scores;

    // define the GA population
    sycl_wrapper<std::vector, chromosome> chromosomes;
    sycl_wrapper<std::vector, chromosome> best_chromosomes;

    // Grid Maps
    // Can be moved inside shared memory for example
    const index3D index_maps;
    const point3D center, minimum, maximum;
    // TODO maps and textures
    sycl_object<fp_type> electro_tex, desolv_tex;
    wrappers_container atom_texs;

    // Random generation
    sycl_random_object random_states;

  public:
    virtual_screen_sycl(const knobs k,
                        std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                        std::shared_ptr<const grid_map>& electro_map,
                        std::shared_ptr<const grid_map>& desolv_map,
                        sycl::queue& queue);

    void operator()(batch& incoming_batch);
  };
} // namespace mudock
