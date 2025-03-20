#pragma once

#include <mudock/format/ob_wrapper.hpp>
#include <mudock/molecule.hpp>
#include <string_view>

namespace mudock {

  class pdbqt {
  public:
    template<class molecule_type>
      requires is_molecule<molecule_type>
    void parse(molecule_type&& molecule, std::string_view molecule_description) const;
  };

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<class molecule_type>
    requires is_molecule<molecule_type>
  void pdbqt::parse(molecule_type&& molecule, std::string_view molecule_description) const {
    // NOTE: we use OpenBabel to parse everything and deal with all the chemistry complexity
    convert(molecule, parse_pdbqt(molecule_description));
  }
} // namespace mudock
