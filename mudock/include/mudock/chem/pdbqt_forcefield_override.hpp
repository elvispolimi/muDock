#pragma once

#include <filesystem>
#include <mudock/chem/autodock_layer.hpp>
#include <mudock/format/format_id.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<class autodock_molecule_type, class source_molecule_type>
    requires derived_from_autodock_layer<autodock_molecule_type> && is_molecule<source_molecule_type>
  void check_source_path_pdbqt_forcefield_override(autodock_molecule_type& target,
                                                   const source_molecule_type& source_molecule) {
    const auto& source_path = source_molecule.properties.get(property_type::SOURCE_PATH);
    if (source_path == "N/A")
      return;

    const auto input_path = std::filesystem::path{source_path};
    if (input_path.has_extension() && parse_supported_format(input_path) == supported_format::PDBQT)
      apply_autodock_forcefield_pdbqt(target, input_path);
  }
} // namespace mudock
