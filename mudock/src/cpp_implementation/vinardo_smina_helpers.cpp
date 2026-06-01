#include <mudock/chem/vinardo_smina_helpers.hpp>

#include <stdexcept>

namespace mudock {

//Here we the helper that are needed to replicate the smina pipeline for scoring, so
// the mobility matrix based on the parsed tree structure
// the num of torsional degrees that are used in the affinity function, in which are used to normalize the result

void mark_not_variable(std::vector<std::uint8_t>& variable,
                       const std::size_t num_atoms,
                       const std::size_t i,
                       const std::size_t j) {
  variable[i * num_atoms + j] = 0;
  variable[j * num_atoms + i] = 0;
}

void apply_smina_fixed_marks(const pdbqt_torsion_branch_node& node,
                             std::vector<std::uint8_t>& variable,
                             const std::size_t num_atoms) {
  //mark as non movable all the atoms inside the same branch
  for (std::size_t i = 0; i < node.atom_indices.size(); ++i) {
    for (std::size_t j = i + 1; j < node.atom_indices.size(); ++j) {
      mark_not_variable(variable, num_atoms, node.atom_indices[i], node.atom_indices[j]);
    }
  }
  //Iterate over the child branches
  for (const auto& child: node.children) {
    mark_not_variable(variable, num_atoms, child.from_atom_index, child.to_atom_index);
    for (const auto atom_index: child.child->atom_indices) {
      //!!! Note here that with the smina semantics we are excluding as relatively movable
      //couples of atoms that relate to from&to, instead with the fragments
      //this were included
      mark_not_variable(variable, num_atoms, child.from_atom_index, atom_index);
      mark_not_variable(variable, num_atoms, child.to_atom_index, atom_index);
    }
    apply_smina_fixed_marks(*child.child, variable, num_atoms);
  }
}

std::vector<std::uint8_t> build_smina_mobility_matrix(const static_molecule& ligand,
                                                      const pdbqt_torsion_tree& tree) {
  //The process to assign mobility is the following:
  //1) We start with a matrix of 1, where all the atoms are considered movable
  //2) We mark as non movable with respect to the same atom
  //3) We mark as non movable all the atoms that are in the same branch
  //4) For each child branch, we mark as non movable:
  //   - the rotor axis pair from-to;
  //   - from with every atom inside the child branch node;
  //   - to with every atom inside the child branch node.
  //
  //This follows smina's torsion-tree semantics: those distances are treated as fixed
  //even if the standard muDock fragment logic would classify some of them as movable.

  const auto num_atoms = static_cast<std::size_t>(ligand.num_atoms());
  std::vector<std::uint8_t> variable(num_atoms * num_atoms, 1);
  //Adjust the diagonal
  for (std::size_t i = 0; i < num_atoms; ++i) {
    variable[i * num_atoms + i] = 0;
  }

  apply_smina_fixed_marks(tree.root, variable, num_atoms);
  return variable;
}

unsigned smina_num_tors(const static_molecule& ligand,
                        std::span<const pdbqt_rotor> rotors,
                        std::span<const vinardo_atom_type> ligand_types) {
  //Here we count torsions as smina does.
  //First we count, for each atom, how many bonded atoms are non-hydrogen(so heavy)
  //Then we consider it a rotatable valid only if both the atoms of the couple differ from hydrogen
  //and they have both at least 2 bonded heavy atoms

  const auto num_atoms = static_cast<std::size_t>(ligand.num_atoms());
  if (ligand_types.size() != num_atoms) {
    throw std::runtime_error("Ligand type count does not match ligand atom count");
  }

  std::vector<unsigned> bonded_heavy_atoms(num_atoms, 0);
  for (const auto& bond: ligand.get_bonds()) {
    const auto source = static_cast<std::size_t>(bond.source);
    const auto dest   = static_cast<std::size_t>(bond.dest);
    if (!is_hydrogen(ligand_types[dest])) {
      ++bonded_heavy_atoms[source];
    }
    if (!is_hydrogen(ligand_types[source])) {
      ++bonded_heavy_atoms[dest];
    }
  }

  unsigned num_tors = 0;
  for (const auto& rotor: rotors) {
    const auto from = rotor.from_atom_index;
    const auto to   = rotor.to_atom_index;
    if (from >= num_atoms || to >= num_atoms) {
      throw std::runtime_error("PDBQT torsion tree rotor index out of ligand bounds");
    }
    if (is_hydrogen(ligand_types[from]) || is_hydrogen(ligand_types[to])) {
      continue;
    }
    if (bonded_heavy_atoms[from] <= 1 || bonded_heavy_atoms[to] <= 1) {
      continue;
    }
    ++num_tors;
  }
  return num_tors;
}

} // namespace mudock
