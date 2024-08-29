#pragma once

#include <mudock/grid.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/fragments.hpp>
#include <vector>

namespace mudock {

  fp_type calc_energy(const std::span<fp_type> ligand_x,
                      const std::span<fp_type> ligand_y,
                      const std::span<fp_type> ligand_z,
                      const std::span<fp_type> ligand_vol,
                      const std::span<fp_type> ligand_solpar,
                      const std::span<fp_type> ligand_charge,
                      const std::span<std::size_t> ligand_num_hbond,
                      const std::span<fp_type> ligand_Rij_hb,
                      const std::span<fp_type> ligand_Rii,
                      const std::span<fp_type> ligand_epsij_hb,
                      const std::span<fp_type> ligand_epsii,
                      const std::span<autodock_ff> ligand_autodock_type,
                      const std::span<const bond> ligand_bond,
                      const std::size_t num_atoms,
                      const std::size_t num_bonds,
                      const fragments<static_containers>& ligand_fragments,
                      const grid_atom_mapper& grid_maps,
                      const grid_map& electro_map,
                      const grid_map& desolv_map);
}
