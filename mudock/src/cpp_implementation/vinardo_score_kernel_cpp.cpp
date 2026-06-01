#include <mudock/cpp_implementation/vinardo_score_kernel_cpp.hpp>
#include <mudock/compute/vinardo_affinity.hpp>
#include <mudock/compute/vinardo_scoring_function.hpp>

#include <cmath>
#include <limits>

namespace mudock {

  void calc_vinardo_energy(const int batch_atoms,
                           const int batch_ligands,
                           const int scores_per_ligand,
                           const fp_type* __restrict__ x_scratch_b,
                           const fp_type* __restrict__ y_scratch_b,
                           const fp_type* __restrict__ z_scratch_b,
                           const fp_type* __restrict__ prot_x_b,
                           const fp_type* __restrict__ prot_y_b,
                           const fp_type* __restrict__ prot_z_b,
                           const int* __restrict__ vinardo_num_tors_b,
                           const int* __restrict__ pl_offsets_b,
                           const int* __restrict__ pl_counts_b,
                           const int* __restrict__ pl_protein_atom_idx_b,
                           const int* __restrict__ pl_ligand_atom_idx_b,
                           const fp_type* __restrict__ pl_radius_sum_b,
                           const std::uint8_t* __restrict__ pl_hydrophobic_possible_b,
                           const std::uint8_t* __restrict__ pl_hbond_possible_b,
                           const int* __restrict__ ll_offsets_b,
                           const int* __restrict__ ll_counts_b,
                           const int* __restrict__ ll_atom_i_idx_b,
                           const int* __restrict__ ll_atom_j_idx_b,
                           const fp_type* __restrict__ ll_radius_sum_b,
                           const std::uint8_t* __restrict__ ll_hydrophobic_possible_b,
                           const std::uint8_t* __restrict__ ll_hbond_possible_b,
                           fp_type* __restrict__ inter_scores_b,
                           fp_type* __restrict__ intra_scores_b,
                           fp_type* __restrict__ scores_b) {

    //Pass through every ligand
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      //Here we are passing through every ligand
      //So once we have fixed the ligand we have that PL and LL are fixed so reusable for every pose
      //What changes is the spatial dispositions of the atoms of each couple and that
      //is going to change the value of the score
      const int pl_offset = pl_offsets_b[ligand_index];
      const int pl_count  = pl_counts_b[ligand_index];

      const int ll_offset = ll_offsets_b[ligand_index];
      const int ll_count  = ll_counts_b[ligand_index];
      //In the multi-pose context we must keep track of the best pose to apply the correction factor
      //During affinity calculation
      int best_pose_index    = 0;
      fp_type best_raw_score = std::numeric_limits<fp_type>::max();

      //Pass through every pose of the ligand
      for (int pose_index{0}; pose_index < scores_per_ligand; ++pose_index) {
        //For instance given the x coordinate in the scratch we have this memory layout
        //For each ligand we have all the poses stored contiguously
        // ligand 0
        //   pose 0: atom 0, atom 1, atom 2,
        //   pose 1: atom 0, atom 1, atom 2,

        // ligand 1
        //   pose 0: atom 0, atom 1, atom 2,
        //   pose 1: atom 0, atom 1, atom 2,
        //batch_ligands = How many ligands in the batch
        //scores_per_ligand = how many poses per ligand
        //batch_atoms = padding for every pose (max between all the ligands)
        //So to obtain the offset for ligand_index and pose_index we have to:
        // off = ligand_index * batch_atoms * scores_per_ligand + pose_index * batch_atoms

        const int coord_offset = ligand_index * batch_atoms * scores_per_ligand + pose_index * batch_atoms;

        const fp_type* ligand_x = x_scratch_b + coord_offset;
        const fp_type* ligand_y = y_scratch_b + coord_offset;
        const fp_type* ligand_z = z_scratch_b + coord_offset;

        //Now we need to compute for every couple of the pose the score
        //Starting from the PL couples
        fp_type pl_score = fp_type{0};
        for (int pair_index{0}; pair_index < pl_count; ++pair_index) {
          const int pair_offset = pl_offset + pair_index;
          const int protein_atom = pl_protein_atom_idx_b[pair_offset];
          const int ligand_atom  = pl_ligand_atom_idx_b[pair_offset];

          const fp_type dx       = prot_x_b[protein_atom] - ligand_x[ligand_atom];
          const fp_type dy       = prot_y_b[protein_atom] - ligand_y[ligand_atom];
          const fp_type dz       = prot_z_b[protein_atom] - ligand_z[ligand_atom];
          const fp_type distance = std::sqrt(dx * dx + dy * dy + dz * dz);

          if (distance >= fp_type{8}) {
            continue;
          }

          const fp_type surface_distance = distance - pl_radius_sum_b[pair_offset];
          pl_score += compute_vinardo_pair_energy(surface_distance,
                                                  pl_hydrophobic_possible_b[pair_offset],
                                                  pl_hbond_possible_b[pair_offset]);
        }

        //The same for the LL couples
        fp_type ll_score = fp_type{0};
        for (int pair_index{0}; pair_index < ll_count; ++pair_index) {
          const int pair_offset = ll_offset + pair_index;
          const int atom_i      = ll_atom_i_idx_b[pair_offset];
          const int atom_j      = ll_atom_j_idx_b[pair_offset];

          const fp_type dx       = ligand_x[atom_i] - ligand_x[atom_j];
          const fp_type dy       = ligand_y[atom_i] - ligand_y[atom_j];
          const fp_type dz       = ligand_z[atom_i] - ligand_z[atom_j];
          const fp_type distance = std::sqrt(dx * dx + dy * dy + dz * dz);

          if (distance >= fp_type{8}) {
            continue;
          }

          const fp_type surface_distance = distance - ll_radius_sum_b[pair_offset];
          ll_score += compute_vinardo_pair_energy(surface_distance,
                                                  ll_hydrophobic_possible_b[pair_offset],
                                                  ll_hbond_possible_b[pair_offset]);
        }

        // With one pose, affinity is computed directly from the protein-ligand score.
        // With multiple poses, Vina/Vinardo uses the intramolecular score of the best pose
        // as reference for the conformation-dependent correction.
        const int score_offset = ligand_index * scores_per_ligand + pose_index;
        inter_scores_b[score_offset] = pl_score;
        intra_scores_b[score_offset] = ll_score;
        //Keep track of the best pose
        const fp_type raw_score = pl_score + ll_score;
        if (raw_score < best_raw_score) {
          best_raw_score  = raw_score;
          best_pose_index = pose_index;
        }

      }

      const int best_score_offset         = ligand_index * scores_per_ligand + best_pose_index;
      const fp_type reference_intra_score = intra_scores_b[best_score_offset];

      for (int pose_index{0}; pose_index < scores_per_ligand; ++pose_index) {
        const int score_offset = ligand_index * scores_per_ligand + pose_index;
        //Score pre affinity(as in the paper we correct with the reference intra score term)
        const fp_type corrected_score = pose_index == best_pose_index
                                            ? inter_scores_b[score_offset]
                                            : inter_scores_b[score_offset] + intra_scores_b[score_offset] -
                                                  reference_intra_score;

        //the memory layout of scores_b is the following:
        //ligand 0: pose 0, pose 1, pose 2,
        // ligand 1: pose 0, pose 1, pose 2,
        // ligand 2: pose 0, pose 1, pose 2,
        //So to write for the given pose of the ligand ligand index we have to:
        // off = ligand_index * scores_per_ligand + pose_index
        scores_b[score_offset] = vinardo_affinity(corrected_score, vinardo_num_tors_b[ligand_index]);
      }
    }
  }

  template<>
  void vinardo_score_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<this->vinardo_region_name>(calc_vinardo_energy,
                                                batch_atoms,
                                                batch_ligands,
                                                scores_per_ligand,
                                                x_scratch_b,
                                                y_scratch_b,
                                                z_scratch_b,
                                                prot_x_b,
                                                prot_y_b,
                                                prot_z_b,
                                                vinardo_num_tors_b,
                                                pl_offsets_b,
                                                pl_counts_b,
                                                pl_protein_atom_idx_b,
                                                pl_ligand_atom_idx_b,
                                                pl_radius_sum_b,
                                                pl_hydrophobic_possible_b,
                                                pl_hbond_possible_b,
                                                ll_offsets_b,
                                                ll_counts_b,
                                                ll_atom_i_idx_b,
                                                ll_atom_j_idx_b,
                                                ll_radius_sum_b,
                                                ll_hydrophobic_possible_b,
                                                ll_hbond_possible_b,
                                                inter_scores_b,
                                                intra_scores_b,
                                                scores_b);
  }
} // namespace mudock
