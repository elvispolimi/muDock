#pragma once

#include "openbabel/mol.h"

#include <cassert>
#include <cmath>
#include <memory>
#include <mudock/chem.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/data.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

namespace mudock {

  //===------------------------------------------------------------------------------------------------------
  // Useful type alias to work with OpenBabel data structures
  //===------------------------------------------------------------------------------------------------------

  // define a type alias to prevent memory leaks
  using ob_mol_wrapper = std::unique_ptr<OpenBabel::OBMol>;

  //===------------------------------------------------------------------------------------------------------
  // Utility functions to parse convert simple data from OpenBabel to our simple data structure
  //===------------------------------------------------------------------------------------------------------

  // the list of functions that we use to parse a molecule
  [[nodiscard]] ob_mol_wrapper parse_pdbqt(const std::string_view description);
  [[nodiscard]] bond_type parse_ob_bond_type(const OpenBabel::OBBond& bond_type);

  //===------------------------------------------------------------------------------------------------------
  // Translate an OpenBabel molecule to our internal format
  //===------------------------------------------------------------------------------------------------------

  template<class molecule_type>
    requires is_molecule<molecule_type>
  void convert(molecule_type&& dest, const ob_mol_wrapper& source) {
    const size_t num_atoms = source->NumAtoms();
    const size_t num_bonds = source->NumBonds();
    // set the molecule geometry
    // NOTE: for static molecules we need to enforce the constraint on the maximum number
    //       ot atoms or bonds by throwing an exception
    // NOTE: this is set to true due to some mismatch between how rdkit counts heavy bonds, and the actual number of bonds
    if constexpr (std::same_as<std::remove_cvref_t<molecule_type>, static_molecule>) {
      if (num_atoms > max_static_atoms() || num_bonds > max_static_bonds()) {
        throw std::runtime_error("Number of atoms or bonds exceeding static storage");
      }
    }
    dest.resize(num_atoms, num_bonds);

    // TODO check charges
    size_t mudock_atom_index{0};
    std::unordered_map<unsigned int, int> index_translator;
    OpenBabel::OBTypeTable ttab;
    ttab.SetFromType("INT");
    ttab.SetToType("XYZ");
    [[maybe_unused]] unsigned long max_atom_index{0};
    for (auto atom_it = source->BeginAtoms(); atom_it < source->EndAtoms(); ++atom_it) {
      const auto atom         = *atom_it;
      const auto atom_id      = atom->GetId();
      const auto atom_element = parse_element_symbol(ttab.Translate(atom->GetType()));
      // Get which atoms are aromatic
      // TODO check the cast between bool and uinfast8_t
      dest.elements(mudock_atom_index)    = atom_element;
      dest.is_aromatic(mudock_atom_index) = atom->IsAromatic();
      dest.x(mudock_atom_index)           = static_cast<fp_type>(atom->GetX());
      dest.y(mudock_atom_index)           = static_cast<fp_type>(atom->GetY());
      dest.z(mudock_atom_index)           = static_cast<fp_type>(atom->GetZ());
      dest.charge(mudock_atom_index)      = atom->GetPartialCharge();
      index_translator.emplace(atom_id, mudock_atom_index);
      ++mudock_atom_index;
      max_atom_index = std::max(max_atom_index, atom_id);
    }
    // Verify that there are no gaps in molecule indexes
    assert((max_atom_index + 1) == num_atoms == mudock_atom_index);

    // fill the bond information
    auto mudock_bond_index = int{0};
    for (auto bond_it = source->BeginBonds(); bond_it < source->EndBonds(); ++bond_it) {
      const auto bond          = *bond_it;
      const int atom_id_source = bond->GetBeginAtomIdx();
      const int atom_id_dest   = bond->GetEndAtomIdx();
      auto& mudock_bond        = dest.bonds(mudock_bond_index);
      mudock_bond.source       = index_translator.at(atom_id_source);
      mudock_bond.dest         = index_translator.at(atom_id_dest);
      mudock_bond.type         = parse_ob_bond_type(*bond);
      mudock_bond.can_rotate   = bond->IsRotor();

      ++mudock_bond_index;
    }

    // store the molecule name
    auto name = source->GetTitle();
    dest.properties.assign(property_type::NAME, name);
  }

} // namespace mudock
