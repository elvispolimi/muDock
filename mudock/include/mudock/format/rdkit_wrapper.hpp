#pragma once

#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <map>
#include <memory>
#include <mudock/chem.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <stdexcept>
#include <string>
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
  [[nodiscard]] rw_mol_wrapper parse_mol2(const std::string_view description);
  [[nodiscard]] rw_mol_wrapper parse_pdb(const std::string_view description);
  [[nodiscard]] bond_type parse_rdkit_bond_type(const RDKit::Bond::BondType bond_type);

  //===------------------------------------------------------------------------------------------------------
  // Translate an RDKit molecule to our internal format
  //===------------------------------------------------------------------------------------------------------

  template<class molecule_type>
    requires is_molecule<molecule_type>
  void convert(molecule_type&& dest, const rw_mol_wrapper& source) {
    // set the molecule geometry
    // NOTE: for static molecules we need to enforce the constraint on the maximum number
    //       ot atoms or bonds by throwing an exception
    // NOTE: this is set to true due to some mismatch between how rdkit counts heavy bonds, and the actual number of bonds
    static constexpr auto only_heavy_bonds = true;
    if constexpr (std::same_as<std::remove_cvref_t<molecule_type>, static_molecule>) {
      if (source->getNumAtoms() > max_static_atoms() ||
          source->getNumBonds(only_heavy_bonds) > max_static_bonds()) {
        throw std::runtime_error("Number of atoms or bonds exceeding static storage");
      }
    }
    dest.resize(source->getNumAtoms(), source->getNumBonds(true));

    // compute the Marsilli-Gasteiger partial charges for the molecule
    computeGasteigerCharges(*source);

    // fill the atom information (we assume a single conformation)
    // NOTE: we need to store the mapping between our atom index and the rdkit one
    std::unordered_map<unsigned int, std::size_t> index_translator;
    assert(source->getNumConformers() == 1);
    const auto conformation = source->getConformer(0);
    auto mudock_atom_index  = std::size_t{0};
    for (const auto& atom: source->atoms()) {
      const auto atom_id      = atom->getIdx();
      const auto [x, y, z]    = conformation.getAtomPos(atom_id);
      const auto atom_element = parse_element_symbol(atom->getSymbol());
      // Get which atoms are aromatic
      // TODO check the cast between bool and uinfast8_t
      assert(atom_element.has_value());
      dest.elements(mudock_atom_index)    = atom_element.value();
      dest.is_aromatic(mudock_atom_index) = atom->getIsAromatic();
      dest.x(mudock_atom_index)           = static_cast<fp_type>(x);
      dest.y(mudock_atom_index)           = static_cast<fp_type>(y);
      dest.z(mudock_atom_index)           = static_cast<fp_type>(z);
      dest.charge(mudock_atom_index)      = fp_type{0};
      std::string charge = atom->getProp<std::string>(RDKit::common_properties::_GasteigerCharge);
      if (charge == "inf")
        // https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/MolSurf.py line 213
        dest.charge(mudock_atom_index) = 0.0;
      else
        dest.charge(mudock_atom_index) = std::stof(charge);
      index_translator.emplace(atom_id, mudock_atom_index);
      mudock_atom_index += std::size_t{1};
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
    auto mudock_bond_index = std::size_t{0};
    for (const auto& rdkit_bond: source->bonds()) {
      const auto atom_id_source = rdkit_bond->getBeginAtomIdx();
      const auto atom_id_dest   = rdkit_bond->getEndAtomIdx();
      auto& bond                = dest.bonds(mudock_bond_index);
      bond.source               = index_translator.at(atom_id_source);
      bond.dest                 = index_translator.at(atom_id_dest);
      bond.type                 = parse_rdkit_bond_type(rdkit_bond->getBondType());

      // check if it can rotate
      for (const auto& rotatable_bond: matched_bonds) {
        assert(rotatable_bond.size() == std::size_t{2}); // make sure that we have a single bond
        const auto atom_id_1 = index_translator.at(rotatable_bond[0].second);
        const auto atom_id_2 = index_translator.at(rotatable_bond[1].second);
        if ((atom_id_source == atom_id_1 && atom_id_dest == atom_id_2) ||
            (atom_id_source == atom_id_2 && atom_id_dest == atom_id_1)) {
          bond.can_rotate = true;
          break;
        }
      }
      ++mudock_bond_index;
    }

    // store the molecule name
    auto name = std::string{"N/A"};
    source->getPropIfPresent<std::string>(RDKit::common_properties::_Name, name);
    dest.properties.assign(property_type::NAME, name);
  }

} // namespace mudock
