#pragma once

#include <mudock/format/rdkit_wrapper.hpp>
#include <mudock/molecule.hpp>
#include <string_view>

namespace mudock {

  class pdb {
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
  void pdb::parse(molecule_type&& molecule, std::string_view molecule_description) const {
    // NOTE: we use RDKit to parse everything and deal with all the chemistry complexity
    convert(molecule, parse_pdb(molecule_description));
  }

} // namespace mudock
