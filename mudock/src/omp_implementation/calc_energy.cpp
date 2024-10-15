#include <cuda.h>
#include <mudock/chem.hpp>
#include <mudock/omp_implementation/calc_energy.hpp>

namespace mudock {
  // TODO should be precomputed some value
  // TODO should be moved into a header file
  fp_type calc_ddd_Mehler_Solmajer_omp(fp_type distance) {
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

  fp_type calc_intra_energy(const fp_type* ligand_x,
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
                            const int ligand_num_nonbonds,
                            const int* __restrict__ ligand_nonbond_a1,
                            const int* __restrict__ ligand_nonbond_a2) {
    fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
    // TODO OPT: precompute value from branch rework_cpp
    // TODO
    return elect_total_eintcal + emap_total_eintcal + dmap_total_eintcal;
  }
} // namespace mudock
