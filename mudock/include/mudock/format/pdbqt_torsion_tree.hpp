#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string_view>
#include <vector>

namespace mudock {

//Helper and struct needed to build a tree representation of the ligand in the PDBQT file.
//It's necessary to build the mobility matrix and to retrieve the num of rotatable bonds.

struct pdbqt_torsion_branch_node;

struct pdbqt_torsion_branch {
  std::size_t from_atom_index;
  std::size_t to_atom_index;
  pdbqt_torsion_branch_node* child;
};

struct pdbqt_torsion_branch_node {
  std::vector<std::size_t> atom_indices;
  std::vector<pdbqt_torsion_branch> children;
};

struct pdbqt_rotor {
  std::size_t from_atom_index;
  std::size_t to_atom_index;
};

struct pdbqt_ligand_info {
  bool valid = false;
  std::vector<std::uint8_t> mobility_matrix;
  std::vector<pdbqt_rotor> rotors;
};

struct pdbqt_torsion_tree {
  pdbqt_torsion_branch_node root;
  std::vector<pdbqt_torsion_branch_node> storage; //Need to return also the storage, otherwise pointer in the branches would be dangling.
  std::vector<pdbqt_rotor> rotors;
};

[[nodiscard]] pdbqt_torsion_tree parse_pdbqt_torsion_tree(const std::filesystem::path& pdbqt_path);
[[nodiscard]] pdbqt_torsion_tree parse_pdbqt_torsion_tree(std::string_view pdbqt_description);

} // namespace mudock
