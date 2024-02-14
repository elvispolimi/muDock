#pragma once

#include "mudock/chem/bond_types.hpp"

#include <GraphMol/RWMol.h>
#include <cassert>
#include <map>
#include <memory>
#include <mudock/chem/elements.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <string_view>
#include <unordered_map>

namespace mudock {

  // define a type alias to prevent memory leaks
  using rw_mol_wrapper = std::unique_ptr<RDKit::RWMol>;

  // the list of functions that we use to parse a molecule
  rw_mol_wrapper parse_mol2(const std::string_view description);

  // utility functions to parse information from RDKit string to our data structure
  bond_type parse_rdkit_bond_type(const RDKit::Bond::BondType bond_type);

  // main function that translates an RDKit molecule to our internal format
  template<class molecule_type>
    requires is_molecule<molecule_type>
  void convert(molecule_type&& dest, const rw_mol_wrapper& source) {
    // set the molecule geometry
    static constexpr auto only_heavy_bonds = false;
    dest.resize(source->getNumAtoms(), source->getNumBonds(only_heavy_bonds));

    // fill the atom information (we assume a single conformation)
    // NOTE: we need to store the mapping between our atom index and the rdkit one
    std::unordered_map<unsigned int, index_type> index_translator;
    assert(source->getNumConformers() == 1);
    const auto conformation = source->getConformer(0);
    auto mudock_atom_index  = index_type{0};
    for (const auto& atom: source->atoms()) {
      const auto atom_id                    = atom->getIdx();
      const auto [x, y, z]                  = conformation.getAtomPos(atom_id);
      dest.elements[mudock_atom_index]      = parse_element_symbol(atom->getSymbol());
      dest.coordinates.x[mudock_atom_index] = static_cast<coordinate_type>(x);
      dest.coordinates.y[mudock_atom_index] = static_cast<coordinate_type>(y);
      dest.coordinates.z[mudock_atom_index] = static_cast<coordinate_type>(z);
      index_translator.emplace(atom_id, mudock_atom_index);
      mudock_atom_index += index_type{1};
    }

    // fill the bond information
    auto mudock_bond_index = index_type{0};
    for (const auto& bond: source->bonds()) {
      dest.bonds[mudock_bond_index].source = index_translator.at(bond->getBeginAtomIdx());
      dest.bonds[mudock_bond_index].dest   = index_translator.at(bond->getEndAtomIdx());
      dest.bonds[mudock_bond_index].type   = parse_rdkit_bond_type(bond->getBondType());
    }
  }

} // namespace mudock
