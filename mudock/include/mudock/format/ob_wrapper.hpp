#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <memory>
#include <mudock/chem.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/chem/elements.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/data.h>
#include <openbabel/generic.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sys/types.h>

namespace mudock {

  //===------------------------------------------------------------------------------------------------------
  // Useful type alias to work with OpenBabel data structures
  //===------------------------------------------------------------------------------------------------------

  using ob_mol_wrapper = std::unique_ptr<OpenBabel::OBMol>;

  enum class supported_format : int { MOL2 = 0, PDBQT, PDB };

  struct format_description {
    supported_format format;
    std::string_view extension;
  };

  static constexpr std::array<format_description, 3> FORMAT_EXTENSIONS = {
      {{supported_format::MOL2, "mol2"}, {supported_format::PDBQT, "pdbqt"}, {supported_format::PDB, "pdb"}}};

  [[nodiscard]] supported_format parse_supported_format(const std::string_view extension);
  [[nodiscard]] inline supported_format parse_supported_format(const std::filesystem::path path) {
    assert(path.has_extension());
    const auto extension = path.extension();
    assert(!extension.empty());
    return parse_supported_format(std::string_view(extension.string().substr(1)));
  }
  [[nodiscard]] std::string_view parse_supported_format(const supported_format format);

  [[nodiscard]] ob_mol_wrapper parser(const std::filesystem::path file_path);
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
  ob_mol_wrapper format_parser<supported_format::PDB>(const std::string_view description);

  template<supported_format format>
  void format_writer(const ob_mol_wrapper& mol, const std::filesystem::path out_path);

  template<>
  void format_writer<supported_format::PDBQT>(const ob_mol_wrapper& mol,
                                              const std::filesystem::path out_path);
  template<>
  void format_writer<supported_format::MOL2>(const ob_mol_wrapper& mol, const std::filesystem::path out_path);
  template<>
  void format_writer<supported_format::PDB>(const ob_mol_wrapper& mol, const std::filesystem::path out_path);

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
      // TODO check the cast between bool and uinfast8_t
      dest.elements(mudock_atom_index) = atom_element;
      // TODO apparently is not working on PDBQT files, manually sets aromaticity
      // Ask Gadio how can I compute aromaticity!!
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
    assert(((max_atom_index + 1) == num_atoms) && (num_atoms == mudock_atom_index));

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

  template<class molecule_type>
    requires is_molecule<molecule_type>
  void convert(const ob_mol_wrapper& dest, molecule_type&& source) {
    // TODO
    throw std::runtime_error("Converted not implemented yet!");
  }

  template<class molecule_type>
    requires is_molecule<molecule_type>
  void parse(molecule_type&& molecule, const std::filesystem::path input_path) {
    convert(molecule, parser(input_path));
  }

  template<supported_format format, class molecule_type>
    requires is_molecule<molecule_type>
  void parse(molecule_type&& molecule, const std::string_view description) {
    convert(molecule, format_parser<format>(description));
  }

  namespace openbabel {
    // The following function relies on the same order of atom loading and file
    template<class molecule_type>
      requires is_molecule<molecule_type>
    void apply_autodock_forcefield(molecule_type&& molecule, const std::filesystem::path input_path) {
      const auto format = parse_supported_format(input_path);
      assert(format == supported_format::PDBQT);

      const auto desc = read_from_stream(std::ifstream(input_path));
      std::stringstream desc_s{desc};

      const std::size_t num_atoms = molecule.num_atoms();

      std::size_t index = 0;
      std::string line;
      while (std::getline(desc_s, line)) {
        if (line.find(pdbqt::PDBQT_ATOM_TOKEN) != std::string::npos ||
            line.find(pdbqt::PDBQT_HETATOM_TOKEN) != std::string::npos) {
          assert(index < num_atoms);
          // FIXMED Really bad, at the moment we rely on OpenBabel PDBQT structure
          // What if PDBQT is standardized..
          if (line.size() < 79)
            line += " ";
          std::string adt_value = line.substr(77, 2);
          assert(adt_value.size() == 2);
          if (adt_value[1] == ' ')
            adt_value.pop_back();
          // FIX ME add check that the order of atoms is the same

          const auto adt                = parse_autodock_type(adt_value);
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
    }

  } // namespace openbabel
} // namespace mudock
