#pragma once

#include <filesystem>
#include <fstream>
#include <mudock/format/supported_format.hpp>
#include <mudock/molecule.hpp>
#include <mudock/utils.hpp>
#include <stdexcept>

namespace mudock {
  template<supported_format format, class molecule_type>
    requires is_molecule<molecule_type>
  [[nodiscard]] molecule_type parser(std::string_view description);

  template<>
  static_molecule parser<supported_format::PDBQT, static_molecule>(std::string_view);
  template<>
  dynamic_molecule parser<supported_format::PDBQT, dynamic_molecule>(std::string_view);
  template<>
  static_molecule parser<supported_format::MOL2, static_molecule>(std::string_view);
  template<>
  dynamic_molecule parser<supported_format::MOL2, dynamic_molecule>(std::string_view);
  template<>
  static_molecule parser<supported_format::PDB, static_molecule>(std::string_view);
  template<>
  dynamic_molecule parser<supported_format::PDB, dynamic_molecule>(std::string_view);
  template<>
  static_molecule parser<supported_format::SDF, static_molecule>(std::string_view);
  template<>
  dynamic_molecule parser<supported_format::SDF, dynamic_molecule>(std::string_view);
  template<>
  static_molecule parser<supported_format::MOL, static_molecule>(std::string_view);
  template<>
  dynamic_molecule parser<supported_format::MOL, dynamic_molecule>(std::string_view);

  template<class molecule_type>
    requires is_molecule<molecule_type>
  molecule_type parser(const std::filesystem::path& file_path) {
    if (!file_path.has_extension()) throw std::runtime_error("Molecular input has no file extension");
    const auto format = parse_supported_format(file_path);
    const auto description = read_from_stream(std::ifstream(file_path, std::ios::binary));
    molecule_type molecule;
    constexpr_switch<0, get_num_supported_format(), 1>(
        [&](const auto index) {
          constexpr auto selected = static_cast<supported_format>(decltype(index)::value);
          if (selected == format) molecule = parser<selected, molecule_type>(description);
        },
        format);
    molecule.properties.assign(property_type::SOURCE_PATH, file_path.string());
    return molecule;
  }
} // namespace mudock
