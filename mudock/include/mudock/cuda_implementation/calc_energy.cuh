#pragma once

#include <mudock/chem/autodock_types.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
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
                                    fp_type* dmap_total_eintcal);
}