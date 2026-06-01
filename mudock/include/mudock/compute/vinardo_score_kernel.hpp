#pragma once

#include <concepts>
#include <cstdint>
#include <memory>
#include <mudock/compute/queue.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct vinardo_score_kernel {
    static constexpr char vinardo_region_name[] = "vinardo_score_kernel";
    vinardo_score_kernel(const int scores_per_ligand_,
                         const int batch_ligands_,
                         const int batch_atoms_,
                         //Ligands positions
                         const fp_type *__restrict__ x_scratch_b_,
                         const fp_type *__restrict__ y_scratch_b_,
                         const fp_type *__restrict__ z_scratch_b_,
                         //Protein position
                         const fp_type *__restrict__ prot_x_b_,
                         const fp_type *__restrict__ prot_y_b_,
                         const fp_type *__restrict__ prot_z_b_,

                         const int *__restrict__ vinardo_num_tors_b_,
                         //PL pairs
                         const int *__restrict__ pl_offsets_b_,
                         const int *__restrict__ pl_counts_b_,
                         const int *__restrict__ pl_protein_atom_idx_b_,
                         const int *__restrict__ pl_ligand_atom_idx_b_,
                         const fp_type *__restrict__ pl_radius_sum_b_,
                         const std::uint8_t *__restrict__ pl_hydrophobic_possible_b_,
                         const std::uint8_t *__restrict__ pl_hbond_possible_b_,
                         //LL pairs
                         const int *__restrict__ ll_offsets_b_,
                         const int *__restrict__ ll_counts_b_,
                         const int *__restrict__ ll_atom_i_idx_b_,
                         const int *__restrict__ ll_atom_j_idx_b_,
                         const fp_type *__restrict__ ll_radius_sum_b_,
                         const std::uint8_t *__restrict__ ll_hydrophobic_possible_b_,
                         const std::uint8_t *__restrict__ ll_hbond_possible_b_,

                         fp_type *__restrict__ inter_scores_b_,
                         fp_type *__restrict__ intra_scores_b_,
                         fp_type *__restrict__ scores_b_,
                         std::shared_ptr<queue_type> q_)
        : scores_per_ligand(scores_per_ligand_),
          batch_ligands(batch_ligands_),
          batch_atoms(batch_atoms_),
          x_scratch_b(x_scratch_b_),
          y_scratch_b(y_scratch_b_),
          z_scratch_b(z_scratch_b_),
          prot_x_b(prot_x_b_),
          prot_y_b(prot_y_b_),
          prot_z_b(prot_z_b_),
          vinardo_num_tors_b(vinardo_num_tors_b_),
          pl_offsets_b(pl_offsets_b_),
          pl_counts_b(pl_counts_b_),
          pl_protein_atom_idx_b(pl_protein_atom_idx_b_),
          pl_ligand_atom_idx_b(pl_ligand_atom_idx_b_),
          pl_radius_sum_b(pl_radius_sum_b_),
          pl_hydrophobic_possible_b(pl_hydrophobic_possible_b_),
          pl_hbond_possible_b(pl_hbond_possible_b_),
          ll_offsets_b(ll_offsets_b_),
          ll_counts_b(ll_counts_b_),
          ll_atom_i_idx_b(ll_atom_i_idx_b_),
          ll_atom_j_idx_b(ll_atom_j_idx_b_),
          ll_radius_sum_b(ll_radius_sum_b_),
          ll_hydrophobic_possible_b(ll_hydrophobic_possible_b_),
          ll_hbond_possible_b(ll_hbond_possible_b_),
          inter_scores_b(inter_scores_b_),
          intra_scores_b(intra_scores_b_),
          scores_b(scores_b_),
          q(q_) {}

    void operator()();

    vinardo_score_kernel(const vinardo_score_kernel &)            = default;
    vinardo_score_kernel(vinardo_score_kernel &&)                 = default;
    vinardo_score_kernel &operator=(const vinardo_score_kernel &) = delete;
    vinardo_score_kernel &operator=(vinardo_score_kernel &&)      = delete;

    ~vinardo_score_kernel() = default;

  private:
    const int scores_per_ligand;
    const int batch_ligands;
    const int batch_atoms;

    const fp_type *__restrict__ x_scratch_b;
    const fp_type *__restrict__ y_scratch_b;
    const fp_type *__restrict__ z_scratch_b;

    const fp_type *__restrict__ prot_x_b;
    const fp_type *__restrict__ prot_y_b;
    const fp_type *__restrict__ prot_z_b;

    const int *__restrict__ vinardo_num_tors_b;

    const int *__restrict__ pl_offsets_b;
    const int *__restrict__ pl_counts_b;
    const int *__restrict__ pl_protein_atom_idx_b;
    const int *__restrict__ pl_ligand_atom_idx_b;
    const fp_type *__restrict__ pl_radius_sum_b;
    const std::uint8_t *__restrict__ pl_hydrophobic_possible_b;
    const std::uint8_t *__restrict__ pl_hbond_possible_b;

    const int *__restrict__ ll_offsets_b;
    const int *__restrict__ ll_counts_b;
    const int *__restrict__ ll_atom_i_idx_b;
    const int *__restrict__ ll_atom_j_idx_b;
    const fp_type *__restrict__ ll_radius_sum_b;
    const std::uint8_t *__restrict__ ll_hydrophobic_possible_b;
    const std::uint8_t *__restrict__ ll_hbond_possible_b;

    fp_type *__restrict__ inter_scores_b;
    fp_type *__restrict__ intra_scores_b;
    fp_type *__restrict__ scores_b;

    std::shared_ptr<queue_type> q;
  };

} // namespace mudock
