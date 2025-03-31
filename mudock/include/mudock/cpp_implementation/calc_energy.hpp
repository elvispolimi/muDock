#pragma once
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  template<cpu_vectorization vect>
  fp_type calc_energy(const fp_type* __restrict__ ligand_x,
                      const fp_type* __restrict__ ligand_y,
                      const fp_type* __restrict__ ligand_z,
                      const fp_type* __restrict__ ligand_vol,
                      const fp_type* __restrict__ ligand_solpar,
                      const fp_type* __restrict__ ligand_charge,
                      const int* __restrict__ map_ligand_offsets,
                      const int num_atoms,
                      const int n_torsions,
                      const int num_nonbond,
                      const int* __restrict__ non_bond_list_a1,
                      const int* __restrict__ non_bond_list_a2,
                      const fp_type* __restrict__ cA_list,
                      const fp_type* __restrict__ cB_list,
                      const int* __restrict__ xB_list,
                      const fp_type* __restrict__ minimum,
                      const fp_type* __restrict__ maximum,
                      const fp_type* __restrict__ center,
                      const int map_index_x,
                      const int map_index_xy,
                      const fp_type* __restrict__ grid_maps,
                      const fp_type* __restrict__ electro_map,
                      const fp_type* __restrict__ desolv_map);
} // namespace mudock
