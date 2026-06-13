#pragma once

#include <array>
#include <assert.h>
#include <filesystem>
#include <string_view>

namespace mudock {
  enum class supported_format : int { MOL2 = 0, PDBQT, PDB, ADTMOL2 };

  struct format_description {
    supported_format format;
    std::string_view extension;
  };

  static constexpr std::array<format_description, 4> FORMAT_EXTENSIONS = {
      {{supported_format::MOL2, "mol2"},
       {supported_format::PDBQT, "pdbqt"},
       {supported_format::PDB, "pdb"},
       {supported_format::ADTMOL2, "adtmol2"}}};

  static constexpr auto get_num_supported_format() { return FORMAT_EXTENSIONS.size(); }

  [[nodiscard]] supported_format parse_supported_format(const std::string_view extension);
  [[nodiscard]] inline supported_format parse_supported_format(const std::filesystem::path path) {
    assert(path.has_extension());
    const auto extension = path.extension();
    assert(!extension.empty());
    return parse_supported_format(std::string_view(extension.string().substr(1)));
  }
  [[nodiscard]] std::string_view parse_supported_format(const supported_format format);
} // namespace mudock
