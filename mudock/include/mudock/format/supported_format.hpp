#pragma once

#include <mudock/format/format_id.hpp>
#include <mudock/format/adt_mol2.hpp>
#include <mudock/format/mol2.hpp>
#include <mudock/format/pdb.hpp>
#include <mudock/format/pdbqt.hpp>

namespace mudock {
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
  struct type_of_format_t<supported_format::ADTMOL2> {
    using type = adt_mol2;
  };

  template<supported_format T>
  using type_of_format = typename type_of_format_t<T>::type;
} // namespace mudock
