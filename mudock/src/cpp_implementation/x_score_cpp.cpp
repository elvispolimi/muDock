#include <cmath>
#include <mudock/chem/x_score_hb.hpp>
#include <mudock/compute/x_score_terms.hpp>
#include <mudock/cpp_implementation/x_score_cpp.hpp>

namespace mudock {

  static constexpr fp_type x_score_dist_cutoff = fp_type{8.0};

  static constexpr int x_hb_hydrophobic = static_cast<int>(x_score_hb::H);

  // Compute the X-Score terms (VDW and HP) for every ligand in the batch and write them into terms_b.
  inline void calc_x_score(const int batch_atoms,
                           const int batch_ligands,
                           const int *__restrict__ num_atoms_b,
                           const fp_type *__restrict__ lig_x_b,
                           const fp_type *__restrict__ lig_y_b,
                           const fp_type *__restrict__ lig_z_b,
                           const fp_type *__restrict__ lig_vdw_b,
                           const int *__restrict__ lig_scorable_b,
                           const int *__restrict__ lig_hb_b,
                           const fp_type *__restrict__ lig_rt_b,
                           const fp_type *__restrict__ lig_hbt_b,
                           const int num_prot_atoms,
                           const fp_type *__restrict__ prot_x_b,
                           const fp_type *__restrict__ prot_y_b,
                           const fp_type *__restrict__ prot_z_b,
                           const fp_type *__restrict__ prot_vdw_b,
                           const int *__restrict__ prot_scorable_b,
                           const int *__restrict__ prot_hb_b,
                           fp_type *__restrict__ terms_b) {
    for (int ligand_index = 0; ligand_index < batch_ligands; ++ligand_index) {
      const int atom_stride = ligand_index * batch_atoms;
      const int num_atoms   = num_atoms_b[ligand_index];

      const fp_type *__restrict__ lig_x    = lig_x_b + atom_stride;
      const fp_type *__restrict__ lig_y    = lig_y_b + atom_stride;
      const fp_type *__restrict__ lig_z    = lig_z_b + atom_stride;
      const fp_type *__restrict__ lig_vdw  = lig_vdw_b + atom_stride;
      const int *__restrict__ lig_scorable = lig_scorable_b + atom_stride;
      const int *__restrict__ lig_hb       = lig_hb_b + atom_stride;

      fp_type vdw_sum = 0;
      fp_type hp_sum  = 0;

      for (int i = 0; i < num_atoms; ++i) {
        if (!lig_scorable[i])
          continue;

        const fp_type lx           = lig_x[i];
        const fp_type ly           = lig_y[i];
        const fp_type lz           = lig_z[i];
        const fp_type lr           = lig_vdw[i];
        const bool lig_hydrophobic = (lig_hb[i] == x_hb_hydrophobic);

        fp_type vdw_asum = 0;
        fp_type hp_asum  = 0;

        for (int j = 0; j < num_prot_atoms; ++j) {
          if (!prot_scorable_b[j])
            continue;

          const fp_type dx = lx - prot_x_b[j];
          const fp_type dy = ly - prot_y_b[j];
          const fp_type dz = lz - prot_z_b[j];
          const fp_type d  = std::sqrt(dx * dx + dy * dy + dz * dz);

          // van der Waals: (d0/d)^8 - 2*(d0/d)^4 within d <= cutoff 
          if (d <= x_score_dist_cutoff) {
            const fp_type d0   = lr + prot_vdw_b[j];
            fp_type tmp1       = d0 / d;
            tmp1               = tmp1 * tmp1 * tmp1 * tmp1;
            const fp_type tmp2 = tmp1 * tmp1;
            vdw_asum += tmp2 - fp_type{2} * tmp1;
          }

          // hydrophobic pair: linear ramp between two hydrophobic atoms, restricted to d < cutoff
          if (lig_hydrophobic && prot_hb_b[j] == x_hb_hydrophobic && d < x_score_dist_cutoff) {
            const fp_type sum_r = lr + prot_vdw_b[j];
            const fp_type d1    = sum_r + fp_type{0.5f};
            const fp_type d2    = sum_r + fp_type{2.2f};
            if (d < d1)
              hp_asum += fp_type{1};
            else if (d < d2)
              hp_asum += (fp_type{1} / (d1 - d2)) * (d - d2);
          }
        }

        // VDW: flip the sign: positive is favorable
        vdw_asum *= fp_type{-1};
        if (vdw_asum >= fp_type{0})
          vdw_sum += vdw_asum;

        if (lig_hydrophobic)
          hp_sum += hp_asum;
      }

      fp_type *__restrict__ terms = terms_b + ligand_index * static_cast<int>(x_term_count);
      terms[x_term_vdw]           = vdw_sum;
      terms[x_term_hp]            = hp_sum;
      terms[x_term_hb]            = lig_hbt_b[ligand_index];
      terms[x_term_rt]            = lig_rt_b[ligand_index];
      // final XScore's affinity prediction
      terms[x_term_pkd] =
          compute_x_score_pkd(terms[x_term_vdw], terms[x_term_hb], terms[x_term_hp], terms[x_term_rt]);
    }
  };

  template<>
  void x_score_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<x_score_kernel::x_region_name>(calc_x_score,
                                                    batch_atoms,
                                                    batch_ligands,
                                                    num_atoms_b,
                                                    lig_x_b,
                                                    lig_y_b,
                                                    lig_z_b,
                                                    lig_vdw_b,
                                                    lig_scorable_b,
                                                    lig_hb_b,
                                                    lig_rt_b,
                                                    lig_hbt_b,
                                                    num_prot_atoms,
                                                    prot_x_b,
                                                    prot_y_b,
                                                    prot_z_b,
                                                    prot_vdw_b,
                                                    prot_scorable_b,
                                                    prot_hb_b,
                                                    terms_b);
  }
} // namespace mudock