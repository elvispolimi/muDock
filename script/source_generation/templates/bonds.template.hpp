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
  enum class bond_type : std::size_t {{% for bond in data %}
    {@ bond.value @} = {@ bond.index @}, // {@ bond.name @}{% endfor %}
  };

  // this is knowledge about the bond (at least that we use)
  struct bond_type_description {
    bond_type value;
    std::string_view name;
  };
  extern const std::array<bond_type_description, {@num_bonds@}> BOND_DICTIONARY;

  // utility functions to work with them
  inline const bond_type_description& get_description(const bond_type b) {
    assert(BOND_DICTIONARY[static_cast<std::size_t>(b)].value == b);
    return BOND_DICTIONARY[static_cast<std::size_t>(b)];
  }

} // namespace mudock
