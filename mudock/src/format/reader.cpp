

#include <algorithm>
#include <cassert>
#include <fstream>
#include <memory>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/log.hpp>
#include <mudock/utils.hpp>
#include <openbabel/atom.h>
#include <openbabel/babelconfig.h>
#include <openbabel/chargemodel.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/obutil.h>
#include <openbabel/plugin.h>
#include <stdexcept>
#include <string_view>

namespace mudock {
  template<supported_format format>
  ob_mol_wrapper ob_parser(const std::string_view description) {
    OpenBabel::obErrorLog.SetOutputLevel(OpenBabel::obMessageLevel::obError); // Silence everything
    OpenBabel::OBConversion conv;
    const std::string ext{parse_supported_format(format)};
    conv.SetInFormat(ext.c_str());

    std::istringstream desc{std::string(description)};
    auto mol = std::make_unique<OpenBabel::OBMol>();
    if (!conv.Read(mol.get(), &desc)) {
      mudock::error("Couldn't open ", ext, " file.");
      throw std::runtime_error("Parser failed, look to logs for details.");
    }
    return mol;
  }

  // ADT_MOL2
  template<class molecule_type>
    requires is_molecule<molecule_type>
  molecule_type parser_impl_adt_mol2(const std::string_view description) {
    molecule_type mol;
    adt_mol2::parse(mol, description);
    return mol;
  }
  template<>
  dynamic_molecule parser<supported_format::ADTMOL2>(const std::string_view description,
                                                     std::function<bool(OpenBabel::OBBond&)>) {
    return parser_impl_adt_mol2<dynamic_molecule>(description);
  };
  template<>
  static_molecule parser<supported_format::ADTMOL2>(const std::string_view description,
                                                    std::function<bool(OpenBabel::OBBond&)>) {
    return parser_impl_adt_mol2<static_molecule>(description);
  };
  template<>
  ob_mol_wrapper parser<supported_format::ADTMOL2>(const std::string_view description,
                                                   std::function<bool(OpenBabel::OBBond&)>) {
    dynamic_molecule mol = parser_impl_adt_mol2<dynamic_molecule>(description);

    ob_mol_wrapper ob_mol = std::make_unique<OpenBabel::OBMol>();
    convert(ob_mol, mol);
    return ob_mol;
  }

  // Parser custom impl
  template<class molecule_type, supported_format format>
    requires is_molecule<molecule_type>
  molecule_type parser_impl(const std::string_view description,
                            std::function<bool(OpenBabel::OBBond&)> rotor_check) {
    molecule_type mol;
    ob_mol_wrapper ob_mol = parser<format, ob_mol_wrapper>(description);
    convert(mol, ob_mol, rotor_check);
    return mol;
  }
  // PDBQT
  template<>
  static_molecule parser<supported_format::PDBQT>(const std::string_view description,
                                                  std::function<bool(OpenBabel::OBBond&)> rotor_check) {
    return parser_impl<static_molecule, supported_format::PDBQT>(description, rotor_check);
  };
  template<>
  dynamic_molecule parser<supported_format::PDBQT>(const std::string_view description,
                                                   std::function<bool(OpenBabel::OBBond&)> rotor_check) {
    return parser_impl<dynamic_molecule, supported_format::PDBQT>(description, rotor_check);
  };
  template<>
  ob_mol_wrapper parser<supported_format::PDBQT>(const std::string_view description,
                                                 std::function<bool(OpenBabel::OBBond&)>) {
    return ob_parser<supported_format::PDBQT>(description);
  }
  // MOL2
  template<>
  static_molecule parser<supported_format::MOL2>(const std::string_view description,
                                                 std::function<bool(OpenBabel::OBBond&)>) {
    return parser_impl<static_molecule, supported_format::MOL2>(description, ob_rotate_check);
  };
  template<>
  dynamic_molecule parser<supported_format::MOL2>(const std::string_view description,
                                                  std::function<bool(OpenBabel::OBBond&)>) {
    return parser_impl<dynamic_molecule, supported_format::MOL2>(description, ob_rotate_check);
  };
  template<>
  ob_mol_wrapper parser<supported_format::MOL2>(const std::string_view description,
                                                std::function<bool(OpenBabel::OBBond&)>) {
    return ob_parser<supported_format::MOL2>(description);
  }
  // PDBQT
  template<>
  static_molecule parser<supported_format::PDB>(const std::string_view description,
                                                std::function<bool(OpenBabel::OBBond&)> rotor_check) {
    return parser_impl<static_molecule, supported_format::PDB>(description, rotor_check);
  };
  template<>
  dynamic_molecule parser<supported_format::PDB>(const std::string_view description,
                                                 std::function<bool(OpenBabel::OBBond&)> rotor_check) {
    return parser_impl<dynamic_molecule, supported_format::PDB>(description, rotor_check);
  };
  template<>
  ob_mol_wrapper parser<supported_format::PDB>(const std::string_view description,
                                               std::function<bool(OpenBabel::OBBond&)>) {
    return ob_parser<supported_format::PDB>(description);
  }
} // namespace mudock
