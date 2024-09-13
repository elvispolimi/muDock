#pragma once

#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/grid.hpp>
#include <mudock/molecule.hpp>
#include <vector>

namespace mudock {
  fp_type calc_energy(const std::span<fp_type> ligand_x,
                      const std::span<fp_type> ligand_y,
                      const std::span<fp_type> ligand_z,
                      const std::span<fp_type> ligand_vol,
                      const std::span<fp_type> ligand_solpar,
                      const std::span<fp_type> ligand_charge,
                      const std::span<int> ligand_num_hbond,
                      const std::span<fp_type> ligand_Rij_hb,
                      const std::span<fp_type> ligand_Rii,
                      const std::span<fp_type> ligand_epsij_hb,
                      const std::span<fp_type> ligand_epsii,
                      const std::span<autodock_ff> ligand_autodock_type,
                      const int num_atoms,
                      const int n_torsions,
                      const std::span<non_bond_parameter> non_bond_list,
                      const grid_atom_mapper& grid_maps,
                      const grid_map& electro_map,
                      const grid_map& desolv_map);
} // namespace mudock
