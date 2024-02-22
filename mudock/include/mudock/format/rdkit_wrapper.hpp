#pragma once

#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <array>
#include <cassert>
#include <concepts>
#include <map>
#include <memory>
#include <mudock/chem.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

namespace mudock {

  //===------------------------------------------------------------------------------------------------------
  // Useful type alias to work with RDKit data structures
  //===------------------------------------------------------------------------------------------------------

  // define a type alias to prevent memory leaks
  using rw_mol_wrapper = std::unique_ptr<RDKit::RWMol>;

  //===------------------------------------------------------------------------------------------------------
  // Utility functions to parse convert simple data from RDKit to our simple data structure
  //===------------------------------------------------------------------------------------------------------

  // the list of functions that we use to parse a molecule
  rw_mol_wrapper parse_mol2(const std::string_view description);
  rw_mol_wrapper parse_pdb(const std::string_view description);
  bond_type parse_rdkit_bond_type(const RDKit::Bond::BondType bond_type);

  //===------------------------------------------------------------------------------------------------------
  // Translate an RDKit molecule to our internal format
  //===------------------------------------------------------------------------------------------------------

  template<class molecule_type>
    requires is_molecule<molecule_type>
  void convert(molecule_type&& dest, const rw_mol_wrapper& source) {
    // set the molecule geometry
    // NOTE: for static molecules we need to enforce the constraint on the maximum number
    //       ot atoms or bonds by throwing an exception
    static constexpr auto only_heavy_bonds = false;
    if constexpr (std::same_as<std::remove_cvref_t<molecule_type>, static_molecule>) {
      if (source->getNumAtoms() > max_static_atoms() ||
          source->getNumBonds(only_heavy_bonds) > max_static_bonds()) {
        throw std::runtime_error("Number of atoms or bonds exceeding static storage");
      }
    }
    dest.resize(source->getNumAtoms(), source->getNumBonds(only_heavy_bonds));

    // fill the atom information (we assume a single conformation)
    // NOTE: we need to store the mapping between our atom index and the rdkit one
    std::unordered_map<unsigned int, index_type> index_translator;
    assert(source->getNumConformers() == 1);
    const auto conformation = source->getConformer(0);
    auto mudock_atom_index  = index_type{0};
    for (const auto& atom: source->atoms()) {
      const auto atom_id      = atom->getIdx();
      const auto [x, y, z]    = conformation.getAtomPos(atom_id);
      const auto atom_element = parse_element_symbol(atom->getSymbol());
      assert(atom_element.has_value());
      dest.elements[mudock_atom_index]      = atom_element.value();
      dest.coordinates.x[mudock_atom_index] = static_cast<coordinate_type>(x);
      dest.coordinates.y[mudock_atom_index] = static_cast<coordinate_type>(y);
      dest.coordinates.z[mudock_atom_index] = static_cast<coordinate_type>(z);
      index_translator.emplace(atom_id, mudock_atom_index);
      mudock_atom_index += index_type{1};
    }

    // define the SMARTS pattern of the rotatable bonds
    static const auto rotatable_bonds_pattern = rw_mol_wrapper{
        RDKit::SmartsToMol("[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])("
                           "[CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]="
                           "[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-,:;!@[!$"
                           "(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])(["
                           "CH3])[CH3])]")};

    // search all the instances of the fragment in the target molecule
    std::vector<RDKit::MatchVectType> matched_bonds;
    static constexpr auto uniquify             = true;
    static constexpr auto recursionPossible    = true;
    static constexpr auto useChirality         = false;
    static constexpr auto useQueryQueryMatches = false;
    static constexpr auto maxMatches           = 10000; // for the proteins
    RDKit::SubstructMatch(*source,
                          *rotatable_bonds_pattern,
                          matched_bonds,
                          uniquify,
                          recursionPossible,
                          useChirality,
                          useQueryQueryMatches,
                          maxMatches);
    assert(matched_bonds.size() < static_cast<std::size_t>(maxMatches) && "Too many rotatable bonds");

    // fill the bond information
    auto mudock_bond_index = index_type{0};
    for (const auto& bond: source->bonds()) {
      const auto atom_id_source            = bond->getBeginAtomIdx();
      const auto atom_id_dest              = bond->getEndAtomIdx();
      dest.bonds[mudock_bond_index].source = index_translator.at(atom_id_source);
      dest.bonds[mudock_bond_index].dest   = index_translator.at(atom_id_dest);
      dest.bonds[mudock_bond_index].type   = parse_rdkit_bond_type(bond->getBondType());

      // check if it can rotate
      for (const auto& rotatable_bond: matched_bonds) {
        assert(rotatable_bond.size() == std::size_t{2}); // make sure that we have a single bond
        const auto atom_id_1 = index_translator.at(rotatable_bond[0].second);
        const auto atom_id_2 = index_translator.at(rotatable_bond[1].second);
        if ((atom_id_source == atom_id_1 && atom_id_dest == atom_id_2) ||
            (atom_id_source == atom_id_2 && atom_id_dest == atom_id_1)) {
          dest.bonds[mudock_bond_index].can_rotate = true;
          break;
        }
      }
      ++mudock_bond_index;
    }

    // store the molecule name
    auto name = std::string{"N/A"};
    source->getPropIfPresent<std::string>("_Name", name);
    dest.properties.assign(property_type::NAME, name);
  }

} // namespace mudock
