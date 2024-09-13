#include <cuda.h>
#include <cuda_runtime.h>
#include <mudock/chem.hpp>
#include <mudock/cuda_implementation/calc_energy.cuh>

namespace mudock {
  // TODO should be precomputed some value
  __device__ fp_type calc_ddd_Mehler_Solmajer_cuda(fp_type distance) {
    const fp_type lambda{0.003627};
    const fp_type epsilon0{78.4};
    const fp_type A{-8.5525};
    const fp_type B = epsilon0 - A;
    const fp_type rk{7.7839};
    const fp_type lambda_B = -lambda * B;

    fp_type epsilon = A + B / (fp_type{1} + rk * expf(lambda_B * distance));

    if (epsilon < std::numeric_limits<fp_type>::epsilon()) {
      epsilon = fp_type{1.0};
    }
    return epsilon;
  }

  // TODO check math functions if they are from CUDA
  __device__ void calc_intra_energy(const fp_type* ligand_x,
                                    const fp_type* ligand_y,
                                    const fp_type* ligand_z,
                                    const fp_type* ligand_vol,
                                    const fp_type* ligand_solpar,
                                    const fp_type* ligand_charge,
                                    const int* ligand_num_hbond,
                                    const fp_type* ligand_Rij_hb,
                                    const fp_type* ligand_Rii,
                                    const fp_type* ligand_epsij_hb,
                                    const fp_type* ligand_epsii,
                                    const int num_atoms,
                                    const int ligand_num_nonbonds,
                                    const int* __restrict__ ligand_nonbond_a1,
                                    const int* __restrict__ ligand_nonbond_a2,
                                    fp_type* elect_total_eintcal,
                                    fp_type* emap_total_eintcal,
                                    fp_type* dmap_total_eintcal) {
    for (int nonbond_list = threadIdx.x; nonbond_list < ligand_num_nonbonds; nonbond_list += blockDim.x) {
      const int& a1 = ligand_nonbond_a1[nonbond_list];
      const int& a2 = ligand_nonbond_a2[nonbond_list];

      const fp_type distance_two = powf(ligand_x[a1] - ligand_x[a2], 2) +
                                   pow(ligand_y[a1] - ligand_y[a2], 2) + powf(ligand_z[a1] - ligand_z[a2], 2);
      const fp_type distance_two_clamp = std::clamp(distance_two, RMIN_ELEC * RMIN_ELEC, distance_two);
      const fp_type distance           = sqrtf(distance_two_clamp);

      //  Calculate  Electrostatic  Energy
      const fp_type r_dielectric = fp_type{1} / (distance * calc_ddd_Mehler_Solmajer_cuda(distance));
      const fp_type e_elec =
          ligand_charge[a1] * ligand_charge[a2] * ELECSCALE * autodock_parameters::coeff_estat * r_dielectric;
      *elect_total_eintcal += e_elec;

      // Calcuare desolv
      const fp_type nb_desolv = (ligand_vol[a2] * (ligand_solpar[a1] + qsolpar * fabsf(ligand_charge[a1])) +
                                 ligand_vol[a1] * (ligand_solpar[a2] + qsolpar * fabsf(ligand_charge[a2])));

      const fp_type e_desolv = autodock_parameters::coeff_desolv *
                               expf(fp_type{-0.5} / (sigma * sigma) * distance_two_clamp) * nb_desolv;
      *dmap_total_eintcal += e_desolv;

      fp_type e_vdW_Hb{0};
      if (distance_two_clamp < nbc2) {
        const int& hbond_i        = ligand_num_hbond[a1];
        const int& hbond_j        = ligand_num_hbond[a2];
        const fp_type& Rij_hb_i   = ligand_Rij_hb[a1];
        const fp_type& Rij_hb_j   = ligand_Rij_hb[a2];
        const fp_type& Rii_i      = ligand_Rii[a1];
        const fp_type& Rii_j      = ligand_Rii[a2];
        const fp_type& epsij_hb_i = ligand_epsij_hb[a1];
        const fp_type& epsij_hb_j = ligand_epsij_hb[a2];
        const fp_type& epsii_i    = ligand_epsii[a1];
        const fp_type& epsii_j    = ligand_epsii[a2];

        int xA = 12;
        int xB = 6;

        fp_type Rij{0}, epsij{0};
        if ((hbond_i == 1 || hbond_i == 2) && hbond_j > 2) {
          Rij   = Rij_hb_j;
          epsij = epsij_hb_j;
          xB    = 10;
        } else if ((hbond_i > 2) && (hbond_j == 1 || hbond_j == 2)) {
          Rij   = Rij_hb_i;
          epsij = epsij_hb_i;
          xB    = 10;
        } else {
          Rij   = (Rii_i + Rii_j) / fp_type{2};
          epsij = sqrtf(epsii_i * epsii_j);
        }
        if (xA != xB) {
          const fp_type tmp = epsij / (xA - xB);
          const fp_type cA  = tmp * powf(Rij, static_cast<fp_type>(xA)) * xB;
          const fp_type cB  = tmp * powf(Rij, static_cast<fp_type>(xB)) * xA;

          const fp_type rA = powf(distance, static_cast<fp_type>(xA));
          const fp_type rB = powf(distance, static_cast<fp_type>(xB));

          e_vdW_Hb = fminf(EINTCLAMP, (cA / rA - cB / rB));

          if (r_smooth > 0) {
            const fp_type rlow  = distance - r_smooth / 2;
            const fp_type rhigh = distance + r_smooth / 2;
            fp_type energy_smooth{100000};
            for (int j = fmaxf(fp_type{0}, fminf(rlow * A_DIV, fp_type{NEINT - 1}));
                 j <= fminf(fp_type{NEINT - 1}, fminf(rhigh * A_DIV, fp_type{NEINT - 1}));
                 ++j)
              energy_smooth = fminf(energy_smooth, e_vdW_Hb);

            e_vdW_Hb = energy_smooth;
          }
        }
      }
      *emap_total_eintcal += e_vdW_Hb;
    }
  }
} // namespace mudock