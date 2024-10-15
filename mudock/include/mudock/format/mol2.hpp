#pragma once

#include <mudock/format/rdkit_wrapper.hpp>
#include <mudock/molecule.hpp>
#include <string_view>

namespace mudock {

  class mol2 {
  public:
    std::string_view::size_type next_molecule_start_index(std::string_view text) const;

    template<class molecule_type>
      requires is_molecule<molecule_type>
    void parse(molecule_type&& molecule, std::string_view molecule_description) const;
  };

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<class molecule_type>
    requires is_molecule<molecule_type>
  void mol2::parse(molecule_type&& molecule, std::string_view molecule_description) const {
    // NOTE: we use RDKit to parse everything and deal with all the chemistry complexity
    convert(molecule, parse_mol2(molecule_description));
  }

} // namespace mudock
