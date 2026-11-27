#pragma once

#include <concepts>
#include <memory>
#include <mudock/chem/drugscore/drugscore_model.hpp>
#include <mudock/chem/drugscore/drugscore_types.hpp>
#include <mudock/compute/queue.hpp>

namespace mudock {

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct ds_score_kernel {
    static constexpr const char ds_region_name[] = "ds_score_kernel";

    ds_score_kernel(const int batch_atoms_,
                    const int batch_ligands_,
                    const int scores_per_ligand_,
                    const fp_type* __restrict__ batch_ligand_x_,
                    const fp_type* __restrict__ batch_ligand_y_,
                    const fp_type* __restrict__ batch_ligand_z_,
                    const ds_type* __restrict__ batch_ligand_ds_types_,
                    const int* batch_ligand_num_atoms_,
                    const fp_type* __restrict__ protein_x_,
                    const fp_type* __restrict__ protein_y_,
                    const fp_type* __restrict__ protein_z_,
                    const ds_type* __restrict__ protein_ds_types_,
                    const int protein_num_atoms_,
                    const fp_type* __restrict__ ds_scores_,
                    fp_type* __restrict__ scores_b_,
                    std::shared_ptr<queue_type> q_)
        : batch_atoms(batch_atoms_),
          batch_ligands(batch_ligands_),
          scores_per_ligand(scores_per_ligand_),
          batch_ligand_x(batch_ligand_x_),
          batch_ligand_y(batch_ligand_y_),
          batch_ligand_z(batch_ligand_z_),
          batch_ligand_ds_types(batch_ligand_ds_types_),
          batch_ligand_num_atoms(batch_ligand_num_atoms_),
          protein_x(protein_x_),
          protein_y(protein_y_),
          protein_z(protein_z_),
          protein_ds_types(protein_ds_types_),
          protein_num_atoms(protein_num_atoms_),
          ds_scores(ds_scores_),
          scores_b(scores_b_),
          q(q_) {}

    void operator()();

    ds_score_kernel(const ds_score_kernel&)            = default;
    ds_score_kernel(ds_score_kernel&&)                 = default;
    ds_score_kernel& operator=(const ds_score_kernel&) = delete;
    ds_score_kernel& operator=(ds_score_kernel&&)      = delete;

    ~ds_score_kernel() = default;

  private:
    const int batch_atoms;
    const int batch_ligands;
    const int scores_per_ligand;
    const fp_type* __restrict__ batch_ligand_x;
    const fp_type* __restrict__ batch_ligand_y;
    const fp_type* __restrict__ batch_ligand_z;
    const ds_type* __restrict__ batch_ligand_ds_types;
    const int* batch_ligand_num_atoms;
    const fp_type* __restrict__ protein_x;
    const fp_type* __restrict__ protein_y;
    const fp_type* __restrict__ protein_z;
    const ds_type* __restrict__ protein_ds_types;
    const int protein_num_atoms;

    const fp_type* __restrict__ ds_scores;
    fp_type* __restrict__ scores_b;
    std::shared_ptr<queue_type> q;
  };

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct ds_score_precomp_kernel {
    static constexpr const char ds_region_name[] = "ds_score_kernel";

    ds_score_precomp_kernel(const int batch_atoms_,
                            const int batch_ligands_,
                            const int scores_per_ligand_,
                            const fp_type* __restrict__ batch_ligand_x_,
                            const fp_type* __restrict__ batch_ligand_y_,
                            const fp_type* __restrict__ batch_ligand_z_,
                            const ds_type* __restrict__ batch_ligand_ds_types_,
                            const int* batch_ligand_num_atoms_,
                            const fp_type* __restrict__ minimum_,
                            const fp_type* __restrict__ maximum_,
                            const int map_index_x_,
                            const int map_index_xy_,
                            const int map_index_xyz_,
                            const fp_type* __restrict__ grid_maps_,
                            fp_type* __restrict__ scores_b_,
                            std::shared_ptr<queue_type> q_)
        : batch_atoms(batch_atoms_),
          batch_ligands(batch_ligands_),
          scores_per_ligand(scores_per_ligand_),
          batch_ligand_x(batch_ligand_x_),
          batch_ligand_y(batch_ligand_y_),
          batch_ligand_z(batch_ligand_z_),
          batch_ligand_ds_types(batch_ligand_ds_types_),
          batch_ligand_num_atoms(batch_ligand_num_atoms_),
          minimum(minimum_),
          maximum(maximum_),
          map_index_x(map_index_x_),
          map_index_xy(map_index_xy_),
          map_index_xyz(map_index_xyz_),
          grid_maps(grid_maps_),
          scores_b(scores_b_),
          q(q_) {}

    void operator()();

    ds_score_precomp_kernel(const ds_score_precomp_kernel&)            = default;
    ds_score_precomp_kernel(ds_score_precomp_kernel&&)                 = default;
    ds_score_precomp_kernel& operator=(const ds_score_precomp_kernel&) = delete;
    ds_score_precomp_kernel& operator=(ds_score_precomp_kernel&&)      = delete;

    ~ds_score_precomp_kernel() = default;

  private:
    const int batch_atoms;
    const int batch_ligands;
    const int scores_per_ligand;
    const fp_type* __restrict__ batch_ligand_x;
    const fp_type* __restrict__ batch_ligand_y;
    const fp_type* __restrict__ batch_ligand_z;
    const ds_type* __restrict__ batch_ligand_ds_types;
    const int* batch_ligand_num_atoms;
    const fp_type* __restrict__ minimum;
    const fp_type* __restrict__ maximum;
    const int map_index_x;
    const int map_index_xy;
    const int map_index_xyz;
    const fp_type* __restrict__ grid_maps;

    fp_type* __restrict__ scores_b;
    std::shared_ptr<queue_type> q;
  };

} // namespace mudock
