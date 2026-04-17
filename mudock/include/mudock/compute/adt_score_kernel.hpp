#pragma once

#include <concepts>
#include <memory>
#include <mudock/compute/queue.hpp>
#include <mudock/type_alias.hpp>

#include <chrono>  
#include <iostream> 

namespace mudock {
  // TODO check maybe the kernel can be fused togheter with main adt score
  // May become an issue to keep separate the TU and the CUDA/etc dependencies
  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct adt_score_kernel {
    static constexpr char adt_region_name[] = "adt_score_kernel";
    adt_score_kernel(const int scores_per_ligand_,
                     const int batch_ligands_,
                     const int batch_atoms_,
                     const int *__restrict__ num_atoms_b_,
                     const int *__restrict__ num_rotamers_b_,
                     const int *__restrict__ num_nonbonds_b_,
                     const fp_type *__restrict__ x_scratch_b_,
                     const fp_type *__restrict__ y_scratch_b_,
                     const fp_type *__restrict__ z_scratch_b_,
                     const fp_type *__restrict__ vols_b_,
                     const fp_type *__restrict__ solpars_b_,
                     const fp_type *__restrict__ charges_b_,
                     const int *__restrict__ map_offsets_b_,
                     const int *__restrict__ nonbond_a1_b_,
                     const int *__restrict__ nonbond_a2_b_,
                     const fp_type *__restrict__ nonbond_cA_b_,
                     const fp_type *__restrict__ nonbond_cB_b_,
                     const int *__restrict__ nonbond_xB_b_,
                     const fp_type *__restrict__ grid_maps_,
                     const fp_type *__restrict__ minimum_,
                     const fp_type *__restrict__ maximum_,
                     const fp_type *__restrict__ center_,
                     const int map_index_x_,
                     const int map_index_xy_,
                     const int map_index_xyz_,
                     fp_type *__restrict__ scores_b_,
                     std::shared_ptr<queue_type> q_)
        : scores_per_ligand(scores_per_ligand_),
          batch_ligands(batch_ligands_),
          batch_atoms(batch_atoms_),
          num_atoms_b(num_atoms_b_),
          num_rotamers_b(num_rotamers_b_),
          num_nonbonds_b(num_nonbonds_b_),
          x_scratch_b(x_scratch_b_),
          y_scratch_b(y_scratch_b_),
          z_scratch_b(z_scratch_b_),
          vols_b(vols_b_),
          solpars_b(solpars_b_),
          charges_b(charges_b_),
          map_offsets_b(map_offsets_b_),
          nonbond_a1_b(nonbond_a1_b_),
          nonbond_a2_b(nonbond_a2_b_),
          nonbond_cA_b(nonbond_cA_b_),
          nonbond_cB_b(nonbond_cB_b_),
          nonbond_xB_b(nonbond_xB_b_),
          grid_maps(grid_maps_),
          minimum(minimum_),
          maximum(maximum_),
          center(center_),
          map_index_x(map_index_x_),
          map_index_xy(map_index_xy_),
          map_index_xyz(map_index_xyz_),
          scores_b(scores_b_),
          q(q_) {}

    void operator()();

    adt_score_kernel(const adt_score_kernel &)            = default;
    adt_score_kernel(adt_score_kernel &&)                 = default;
    adt_score_kernel &operator=(const adt_score_kernel &) = delete;
    adt_score_kernel &operator=(adt_score_kernel &&)      = delete;

    ~adt_score_kernel() {
        std::cout << "\n[PROFILAZIONE CUSTOM] Tempo TOTALE dentro lo Score Kernel: " 
                  << total_kernel_time << " secondi su " 
                  << total_calls << " chiamate." << std::endl;
    }

  private:
    const int scores_per_ligand;
    const int batch_ligands;
    const int batch_atoms;
    const int *__restrict__ num_atoms_b;
    const int *__restrict__ num_rotamers_b;
    const int *__restrict__ num_nonbonds_b;
    const fp_type *__restrict__ x_scratch_b;
    const fp_type *__restrict__ y_scratch_b;
    const fp_type *__restrict__ z_scratch_b;
    const fp_type *__restrict__ vols_b;
    const fp_type *__restrict__ solpars_b;
    const fp_type *__restrict__ charges_b;
    const int *__restrict__ map_offsets_b;
    const int *__restrict__ nonbond_a1_b;
    const int *__restrict__ nonbond_a2_b;
    const fp_type *__restrict__ nonbond_cA_b;
    const fp_type *__restrict__ nonbond_cB_b;
    const int *__restrict__ nonbond_xB_b;
    const fp_type *__restrict__ grid_maps;
    const fp_type *__restrict__ minimum;
    const fp_type *__restrict__ maximum;
    const fp_type *__restrict__ center;
    const int map_index_x;
    const int map_index_xy;
    const int map_index_xyz;
    fp_type *__restrict__ scores_b;
    std::shared_ptr<queue_type> q;

    //per misurare il tempo
    double total_kernel_time{0.0};
    long long total_calls{0};
  };

} // namespace mudock
