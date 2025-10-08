#pragma once

#include <array>
#include <assert.h>
#include <filesystem>
#include <mudock/format/mol2.hpp>
#include <mudock/format/mol2x.hpp>
#include <mudock/format/pdb.hpp>
#include <mudock/format/pdbqt.hpp>
#include <string>

namespace mudock {
  enum class supported_format : int { MOL2 = 0, PDBQT, PDB, MOL2X };

  // Trait: Kind -> Type
  template<supported_format>
  struct type_of_format_t; // primary (no def)

  template<>
  struct type_of_format_t<supported_format::MOL2> {
    using type = mol2;
  };
  template<>
  struct type_of_format_t<supported_format::PDBQT> {
    using type = pdbqt;
  };
  template<>
  struct type_of_format_t<supported_format::PDB> {
    using type = pdb;
  };
  template<>
  struct type_of_format_t<supported_format::MOL2X> {
    using type = mol2x;
  };

  template<supported_format T>
  using type_of_format = typename type_of_format_t<T>::type;

  struct format_description {
    supported_format format;
    std::string_view extension;
  };

  static constexpr std::array<format_description, 4> FORMAT_EXTENSIONS = {
      {{supported_format::MOL2, "mol2"},
       {supported_format::PDBQT, "pdbqt"},
       {supported_format::PDB, "pdb"},
       {supported_format::MOL2X, "mol2x"}}};

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
