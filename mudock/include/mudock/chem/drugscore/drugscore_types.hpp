#pragma once

#include <array>
#include <cassert>
#include <string_view>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/drugscore_data.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  enum class ds_type : unsigned int {
    Br    = 0,
    C_2   = 1,
    C_3   = 2,
    C_ar  = 3,
    C_cat = 4,
    N_2   = 5,
    N_4   = 6,
    N_am  = 7,
    N_ar  = 8,
    N_pl3 = 9,
    O_2   = 10,
    O_3   = 11,
    O_co2 = 12,
    S_3   = 13,
    C_1   = 14,
    N_3   = 15,
    Cl    = 16,
    F     = 17,
    I     = 18,
    N_1   = 19,
    P_3   = 20,
    S_2   = 21,
    S_o2  = 22,
    S_o   = 23,
    METAL = 24,
    ND    = 25,
  };

  // this is knowledge that we have about all the elements (that we need at least)
  struct ds_type_description {
    ds_type value;
    std::string_view symbol;
  };

  static constexpr unsigned int ds2index(const ds_type e) { return static_cast<unsigned int>(e); }
  static constexpr unsigned int num_ds_types() { return 26; }
  extern const std::array<ds_type_description, num_ds_types()> DS_TYPE_DICTIONARY;

  static constexpr unsigned int num_valued_ds_types() { return 25; }

  // utility functions to work with them
  inline const ds_type_description& get_description(const ds_type e) {
    assert(DS_TYPE_DICTIONARY[static_cast<unsigned int>(e)].value == e);
    return DS_TYPE_DICTIONARY[static_cast<unsigned int>(e)];
  }

  ds_type parse_ds_type_symbol(const std::string_view symbol);

} // namespace mudock
