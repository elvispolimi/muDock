#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <optional>
#include <string_view>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/periodic_table.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  // this is the list of all the the types of intramolecular bond that we know about
  enum class bond_type : int {
    SINGLE   = 0, // Single
    DOUBLE   = 1, // Double
    TRIPLE   = 2, // Triple
    AMIDE    = 3, // Amide
    AROMATIC = 4, // Aromatic
  };

  // this is knowledge about the bond (at least that we use)
  struct bond_type_description {
    bond_type value;
    std::string_view name;
  };
  extern const std::array<bond_type_description, 5> BOND_DICTIONARY;

  // utility functions to work with them
  inline const bond_type_description& get_description(const bond_type b) {
    assert(BOND_DICTIONARY[static_cast<int>(b)].value == b);
    return BOND_DICTIONARY[static_cast<int>(b)];
  }

} // namespace mudock
