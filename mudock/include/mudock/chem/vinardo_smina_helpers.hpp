#pragma once

#include <cstdint>
#include <mudock/chem/vinardo_type.hpp>
#include <mudock/format/pdbqt_torsion_tree.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/containers.hpp>
#include <span>
#include <vector>

namespace mudock {

[[nodiscard]] std::vector<std::uint8_t> build_smina_mobility_matrix(const static_molecule& ligand,const pdbqt_torsion_tree& tree);

[[nodiscard]] unsigned smina_num_tors(const static_molecule& ligand,
                                      std::span<const pdbqt_rotor> rotors,
                                      std::span<const vinardo_atom_type> ligand_types);

} // namespace mudock
