#pragma once

#include "mudock/chem/autodock_ligand.hpp"
#include "mudock/chem/autodock_protein.hpp"

#include <cassert>
#include <cmath>
#include <concepts>
#include <filesystem>
#include <memory>
#include <mudock/chem.hpp>
#include <mudock/chem/autodock_babel_types.hpp>
#include <mudock/chem/autodock_molecule.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/chem/elements.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/data.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/oberror.h>
#include <sys/types.h>
#include <utility>

namespace mudock {

  //===------------------------------------------------------------------------------------------------------
  // Useful type alias to work with OpenBabel data structures
  //===------------------------------------------------------------------------------------------------------

  using ob_mol_wrapper = std::unique_ptr<OpenBabel::OBMol>;

  void writer(const ob_mol_wrapper& mol, const std::filesystem::path out_path);

  //===------------------------------------------------------------------------------------------------------
  // Utility functions to parse convert simple data from OpenBabel to our simple data structure
  //===------------------------------------------------------------------------------------------------------

  // the list of functions that we use to parse a molecule
  template<supported_format format>
  [[nodiscard]] ob_mol_wrapper format_parser(const std::string_view description);

  template<>
  ob_mol_wrapper format_parser<supported_format::PDBQT>(const std::string_view description);
  template<>
  ob_mol_wrapper format_parser<supported_format::MOL2>(const std::string_view description);
  template<>
  ob_mol_wrapper format_parser<supported_format::MOL2X>(const std::string_view description);
  template<>
  ob_mol_wrapper format_parser<supported_format::PDB>(const std::string_view description);

  template<supported_format format, class molecule_type>
    requires derived_from_molecule<molecule_type>
  void format_parser(molecule_type& mol, const std::string_view description);
  template<>
  void format_parser<supported_format::MOL2X, autodock_ligand>(autodock_ligand& mol,
                                                               const std::string_view description);
  template<>
  void format_parser<supported_format::MOL2X, autodock_protein>(autodock_protein& mol,
                                                                const std::string_view description);

  template<supported_format format>
  void format_writer(const ob_mol_wrapper& mol, const std::filesystem::path out_path);
  template<>
  void format_writer<supported_format::PDBQT>(const ob_mol_wrapper& mol,
                                              const std::filesystem::path out_path);
  template<>
  void format_writer<supported_format::MOL2>(const ob_mol_wrapper& mol, const std::filesystem::path out_path);
  template<>
  void format_writer<supported_format::MOL2X>(const ob_mol_wrapper& mol,
                                              const std::filesystem::path out_path);
  template<>
  void format_writer<supported_format::PDB>(const ob_mol_wrapper& mol, const std::filesystem::path out_path);

  template<supported_format format>
  void format_writer(const ob_mol_wrapper& mol, std::ofstream& ofs);
  template<>
  void format_writer<supported_format::MOL2X>(const ob_mol_wrapper& mol, std::ofstream& ofs);
  template<>
  void format_writer<supported_format::PDB>(const ob_mol_wrapper& mol, std::ofstream& ofs);
  template<>
  void format_writer<supported_format::MOL2>(const ob_mol_wrapper& mol, std::ofstream& ofs);
  template<>
  void format_writer<supported_format::PDBQT>(const ob_mol_wrapper& mol, std::ofstream& ofs);

  [[nodiscard]] bond_type parse_ob_bond_type(const OpenBabel::OBBond& bond_type);

  //===------------------------------------------------------------------------------------------------------
  // Translate an OpenBabel molecule to our internal format
  //===------------------------------------------------------------------------------------------------------
  template<typename F>
  concept is_rotate_check = requires(F f, OpenBabel::OBBond& bond) {
    { f(bond) } -> std::same_as<bool>;
  };

  template<auto rotor_check, class molecule_type>
    requires derived_from_molecule<molecule_type> && is_rotate_check<decltype(rotor_check)>
  void convert(molecule_type& dest, const ob_mol_wrapper& source) {
    const size_t num_atoms = source->NumAtoms();
    const size_t num_bonds = source->NumBonds();
    // set the molecule geometry
    // NOTE: for static molecules we need to enforce the constraint on the maximum number
    //       ot atoms or bonds by throwing an exception
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
      const auto atom_type    = atom->GetType();
      const auto atom_element = parse_element_symbol(ttab.Translate(atom_type));
      // Get which atoms are aromatic
      dest.elements(mudock_atom_index) = atom_element;
      dest.x(mudock_atom_index)        = static_cast<fp_type>(atom->GetX());
      dest.y(mudock_atom_index)        = static_cast<fp_type>(atom->GetY());
      dest.z(mudock_atom_index)        = static_cast<fp_type>(atom->GetZ());
      if constexpr (std::derived_from<molecule_type, autodock_static_molecule> ||
                    std::derived_from<molecule_type, autodock_dynamic_molecule>) {
        dest.charge(mudock_atom_index)      = atom->GetPartialCharge();
        dest.is_aromatic(mudock_atom_index) = atom->IsAromatic();
      }
      index_translator.emplace(atom_id, mudock_atom_index);
      ++mudock_atom_index;
      max_atom_index = std::max(max_atom_index, atom_id);
    }
    // Verify that there are no gaps in molecule indexes
    assert(((max_atom_index + 1) == num_atoms) && (num_atoms == mudock_atom_index));

    // fill the bond information
    auto mudock_bond_index = int{0};
    for (auto bond_it = source->BeginBonds(); bond_it < source->EndBonds(); ++bond_it) {
      const auto bond          = *bond_it;
      const int atom_id_source = bond->GetBeginAtomIdx() - 1;
      const int atom_id_dest   = bond->GetEndAtomIdx() - 1;
      auto& mudock_bond        = dest.bonds(mudock_bond_index);
      mudock_bond.source       = index_translator.at(atom_id_source);
      mudock_bond.dest         = index_translator.at(atom_id_dest);
      mudock_bond.type         = parse_ob_bond_type(*bond);
      mudock_bond.can_rotate   = rotor_check(*bond);

      ++mudock_bond_index;
    }

    // store the molecule name
    auto name = source->GetTitle();
    dest.properties.assign(property_type::NAME, name);
  }

  [[nodiscard]] bool rotate_check(::OpenBabel::OBBond&);
  [[nodiscard]] bool pdbqt_rotate_check(::OpenBabel::OBBond&);

  [[nodiscard]] autodock_ff convert_ob_elements(OpenBabel::OBAtom&);

  template<class molecule_type>
    requires is_molecule<molecule_type>
  void apply_autodock_forcefield_ob(molecule_type&& molecule, const ob_mol_wrapper& ob_mol) {
    int index = 0;
    for (auto atom_it = ob_mol->BeginAtoms(); atom_it < ob_mol->EndAtoms(); ++atom_it) {
      const auto atom = *atom_it;

      autodock_ff adt               = convert_ob_elements(*atom);
      molecule.autodock_type(index) = adt;
      const auto& ff_entry          = get_description(adt);
      molecule.autodock_type(index) = ff_entry.value;
      molecule.Rii(index)           = ff_entry.Rii;
      molecule.epsii(index)         = ff_entry.epsii * autodock_parameters::coeff_vdW;
      molecule.vol(index)           = ff_entry.vol;
      molecule.solpar(index)        = ff_entry.solpar;
      molecule.Rij_hb(index)        = ff_entry.Rij_hb;
      molecule.epsij_hb(index)      = ff_entry.epsij_hb * autodock_parameters::coeff_hbond;
      molecule.num_hbond(index)     = ff_entry.hbond;
      ++index;
    }
  }

  template<supported_format format, class molecule_type, auto rotor_check = rotate_check>
    requires derived_from_molecule<molecule_type> && is_rotate_check<decltype(rotor_check)>
  void parse(molecule_type& molecule, const std::string_view description) {
    if constexpr (format == supported_format::MOL2X) {
      format_parser<format, molecule_type>(molecule, description);
    } else {
      const auto ob_mol = format_parser<format>(description);
      convert<rotor_check>(molecule, ob_mol);
    }
  }

  template<class molecule_type, auto rotor_check = rotate_check>
    requires derived_from_molecule<molecule_type> && is_rotate_check<decltype(rotor_check)>
  void parser(molecule_type& molecule, const std::filesystem::path file_path) {
    const auto in_format = parse_supported_format(file_path);

    const auto description = read_from_stream(std::ifstream(file_path));
    constexpr_switch<0, get_num_supported_format(), 1>(
        [&](const auto format_index) {
          const auto format = static_cast<supported_format>(format_index());
          parse<format, molecule_type, rotor_check>(molecule, description);
        },
        in_format);
  }

  template<class molecule_type, auto rotor_check = rotate_check>
    requires derived_from_molecule<molecule_type> && is_rotate_check<decltype(rotor_check)>
  void parse(molecule_type& molecule, const std::filesystem::path input_path) {
    parser<molecule_type, rotor_check>(molecule, input_path);
  }
} // namespace mudock
