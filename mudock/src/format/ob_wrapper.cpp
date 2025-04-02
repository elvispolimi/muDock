#include <algorithm>
#include <cassert>
#include <fstream>
#include <memory>
#include <mudock/format/ob_wrapper.hpp>
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

  supported_format parse_supported_format(const std::string_view extension) {
    const auto element_it = std::find_if(std::begin(FORMAT_EXTENSIONS),
                                         std::end(FORMAT_EXTENSIONS),
                                         [&extension](const auto& e) { return e.extension == extension; });
    if (element_it != std::end(FORMAT_EXTENSIONS))
      return element_it->format;
    else
      throw std::runtime_error("Missing extension");
  }

  std::string_view parse_supported_format(const supported_format format) {
    const auto element_it = std::find_if(std::begin(FORMAT_EXTENSIONS),
                                         std::end(FORMAT_EXTENSIONS),
                                         [&format](const auto& e) { return e.format == format; });
    if (element_it != std::end(FORMAT_EXTENSIONS))
      return element_it->extension;
    else
      throw std::runtime_error("Missing format");
  }

  ob_mol_wrapper parser(const std::filesystem::path file_path) {
    assert(file_path.has_extension());
    const auto extension = file_path.extension();
    assert(!extension.empty());
    const auto format = parse_supported_format(extension.string().substr(1));

    ob_mol_wrapper mol;
    const auto description = read_from_stream(std::ifstream(file_path));
    constexpr_switch<0, FORMAT_EXTENSIONS.size(), 1>(
        [&](const auto format_index) {
          mol = format_parser<static_cast<supported_format>(format_index())>(description);
        },
        format);
    return mol;
  }

  void writer(const ob_mol_wrapper& mol, const std::filesystem::path out_path) {
    assert(out_path.has_extension());
    const auto extension = out_path.extension();
    assert(!extension.empty());
    const auto format = parse_supported_format(extension.string().substr(1));

    constexpr_switch<0, FORMAT_EXTENSIONS.size(), 1>(
        [&](const auto format_index) {
          format_writer<static_cast<supported_format>(format_index())>(mol, out_path);
        },
        format);
  }

  template<supported_format format>
  ob_mol_wrapper ob_parser(const std::string_view description) {
    OpenBabel::OBConversion conv;
    const std::string ext{parse_supported_format(format)};
    conv.SetInFormat(ext.c_str());

    std::istringstream desc{std::string(description)};
    auto mol = std::make_unique<OpenBabel::OBMol>();
    if (!conv.Read(mol.get(), &desc)) {
      mudock::error(std::format("Couldn't open {} file", ext));
      throw std::runtime_error(std::format("{} Parser failed, look to logs for details", ext));
    }
    return mol;
  }

  template<>
  ob_mol_wrapper format_parser<supported_format::PDBQT>(const std::string_view description) {
    auto mol = ob_parser<supported_format::PDBQT>(description);

    return mol;
  }

  template<>
  ob_mol_wrapper format_parser<supported_format::PDB>(const std::string_view description) {
    auto mol = ob_parser<supported_format::PDB>(description);

    return mol;
  }

  template<>
  ob_mol_wrapper format_parser<supported_format::MOL2>(const std::string_view description) {
    auto mol = ob_parser<supported_format::MOL2>(description);
    return mol;
  }

  void format_writer(const ob_mol_wrapper& mol, const std::filesystem::path out_path) {
    assert(out_path.has_extension());
    const auto extension = out_path.extension();
    assert(!extension.empty());
    const auto format = parse_supported_format(extension.string().substr(1));

    constexpr_switch<0, FORMAT_EXTENSIONS.size(), 1>(
        [&](const auto format_index) {
          format_writer<static_cast<supported_format>(format_index())>(mol, out_path);
        },
        format);
  }

  template<supported_format format>
  void ob_writer(const ob_mol_wrapper& mol, const std::filesystem::path out_path) {
    OpenBabel::OBConversion conv;
    const std::string ext{parse_supported_format(format)};
    conv.SetOutFormat(ext.c_str());

    std::ofstream ofs(out_path);
    if (!ofs) {
      throw std::runtime_error("Error: Cannot open file for writing!");
    }

    if (!conv.Write(mol.get(), &ofs)) {
      throw std::runtime_error("Error: Failed to write molecule!");
    }
  }

  template<>
  void format_writer<supported_format::PDBQT>(const ob_mol_wrapper& mol,
                                              const std::filesystem::path out_path) {
    // Add charges using the Gasteiger method
    OpenBabel::OBChargeModel* chargeModel = OpenBabel::OBChargeModel::FindType("gasteiger");
    if (!chargeModel) {
      mudock::error("Error: Unable to find charge model.");
      throw std::runtime_error("Error in OpenBabel transformations");
    }
    chargeModel->ComputeCharges(*mol.get());

    mol.get()->ConnectTheDots();
    mol.get()->PerceiveBondOrders();

    ob_writer<supported_format::PDBQT>(mol, out_path);
  }
  template<>
  void format_writer<supported_format::MOL2>(const ob_mol_wrapper& mol,
                                             const std::filesystem::path out_path) {
    ob_writer<supported_format::MOL2>(mol, out_path);
  }
  template<>
  void format_writer<supported_format::PDB>(const ob_mol_wrapper& mol, const std::filesystem::path out_path) {
    ob_writer<supported_format::PDB>(mol, out_path);
  }

  bond_type parse_ob_bond_type(const OpenBabel::OBBond& bond) {
    if (bond.IsAromatic())
      return bond_type::AROMATIC;
    else if (const_cast<OpenBabel::OBBond&>(bond).IsAmide())
      return bond_type::AMIDE;
    else
      switch (bond.GetBondOrder()) {
        case 1: return bond_type::SINGLE;
        case 2: return bond_type::DOUBLE;
        case 3: return bond_type::TRIPLE;
        default: throw std::runtime_error("Unsopported OpenBabel bond type");
      }
  }

} // namespace mudock
