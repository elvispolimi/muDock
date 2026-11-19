#pragma once

#include <fstream>
#include <mudock/format/mol2.hpp>
#include <mudock/format/supported_format.hpp>
#include <stdexcept>

namespace mudock {
  template<supported_format format, class molecule_type>
    requires is_molecule<molecule_type>
  void writer(const molecule_type& molecule, std::ofstream& output);

  template<>
  inline void writer<supported_format::MOL2, static_molecule>(const static_molecule& molecule,
                                                               std::ofstream& output) {
    mol2::print(molecule, output);
  }
  template<>
  inline void writer<supported_format::MOL2, dynamic_molecule>(const dynamic_molecule& molecule,
                                                                std::ofstream& output) {
    mol2::print(molecule, output);
  }
} // namespace mudock
