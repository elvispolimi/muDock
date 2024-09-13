#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <mudock/chem/elements.hpp>
#include <mudock/type_alias.hpp>
#include <optional>
#include <string_view>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/autodock_babel_types.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  // this is the list of all the known atoms
  enum class autodock_ff: int {{% for element in data %}
    {@ element.value @} = {@ element.index @}, // {@ element.name @}{% endfor %}
  };

  // this is knowledge that we have about all the elements (that we need at least)
  struct autodock_ff_description {
    autodock_ff value;
    std::string_view name;
    fp_type Rii       = 0;
    fp_type epsii     = 0;
    fp_type vol       = 0;
    fp_type solpar    = 0;
    fp_type Rij_hb    = 0;
    fp_type epsij_hb  = 0;
    int hbond = 0;
  };
  extern const std::array < autodock_ff_description, { @num_elements @ }
  > AUTODOCK_FF_DICTIONARY;

  // utility functions to work with them
  inline const autodock_ff_description& get_description(const autodock_ff a) {
    assert(AUTODOCK_FF_DICTIONARY[static_cast<int>(a)].value == a);
    return AUTODOCK_FF_DICTIONARY[static_cast<int>(a)];
  }
} // namespace mudock
