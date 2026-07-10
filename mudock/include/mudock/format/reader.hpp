#pragma once

#include <cassert>
#include <cmath>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/format/supported_format.hpp>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/data.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/oberror.h>
#include <sys/types.h>

namespace mudock {
  template<supported_format format, class molecule>
  [[nodiscard]] molecule parser(const std::string_view description,
                                std::function<bool(OpenBabel::OBBond &)> f = ob_rotate_check);

  template<>
  ob_mol_wrapper parser<supported_format::PDBQT>(const std::string_view description,
                                                 std::function<bool(OpenBabel::OBBond &)>);
  template<>
  ob_mol_wrapper parser<supported_format::MOL2>(const std::string_view description,
                                                std::function<bool(OpenBabel::OBBond &)>);
  template<>
  ob_mol_wrapper parser<supported_format::ADTMOL2>(const std::string_view description,
                                                   std::function<bool(OpenBabel::OBBond &)>);
  template<>
  ob_mol_wrapper parser<supported_format::PDB>(const std::string_view description,
                                               std::function<bool(OpenBabel::OBBond &)>);
  template<>
  static_molecule parser<supported_format::PDBQT>(const std::string_view description,
                                                  std::function<bool(OpenBabel::OBBond &)>);
  template<>
  static_molecule parser<supported_format::MOL2>(const std::string_view description,
                                                 std::function<bool(OpenBabel::OBBond &)>);
  template<>
  static_molecule parser<supported_format::ADTMOL2>(const std::string_view description,
                                                    std::function<bool(OpenBabel::OBBond &)>);
  template<>
  static_molecule parser<supported_format::PDB>(const std::string_view description,
                                                std::function<bool(OpenBabel::OBBond &)>);
  template<>
  dynamic_molecule parser<supported_format::PDBQT>(const std::string_view description,
                                                   std::function<bool(OpenBabel::OBBond &)>);
  template<>
  dynamic_molecule parser<supported_format::MOL2>(const std::string_view description,
                                                  std::function<bool(OpenBabel::OBBond &)>);
  template<>
  dynamic_molecule parser<supported_format::ADTMOL2>(const std::string_view description,
                                                     std::function<bool(OpenBabel::OBBond &)>);
  template<>
  dynamic_molecule parser<supported_format::PDB>(const std::string_view description,
                                                 std::function<bool(OpenBabel::OBBond &)>);

  template<class molecule_type>
    requires is_molecule<molecule_type>
  molecule_type parser(const std::filesystem::path file_path,
                       std::function<bool(OpenBabel::OBBond &)> check_rotor_bond = ob_rotate_check) {
    const auto in_format = parse_supported_format(file_path);
    molecule_type mol;
    const auto description = read_from_stream(std::ifstream(file_path));
    constexpr_switch<0, get_num_supported_format(), 1>(
        [&](const auto format_index) {
          const auto format = static_cast<supported_format>(format_index());
          mol               = parser<format, molecule_type>(description, check_rotor_bond);
        },
        in_format);
    mol.properties.assign(property_type::SOURCE_PATH, file_path.string());
    return mol;
    // TODO
    // mudock::error("The provided path " + file_path.string() + " extension is not yet supported");
    // throw std::runtime_error("Error in parsing");
  }
} // namespace mudock
