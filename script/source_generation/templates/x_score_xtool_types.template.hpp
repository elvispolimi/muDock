#pragma once

#include <array>
#include <cassert>
#include <mudock/type_alias.hpp>
#include <string_view>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/x_score_xtool_types.json and subsequently modified manually
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  // List of all known atoms for X-TOOL force field
  enum class xtool_ff : int {
{% for element in data %}
    {@ element.name @} = {@ element.index @}, // {@ element.name @}
{% endfor %}
  };

  // Knowledge about the X-TOOL force field parameters
  struct xtool_ff_description {
    xtool_ff value;
    std::string_view name;
    fp_type atomic_weight;
    fp_type vdw_radius;
    fp_type vdw_potential;
    fp_type par_charge;
    std::string_view hbond;
  };

  extern const std::array<xtool_ff_description, {@ num_elements @}> XTOOL_FF_DICTIONARY;

  // Utility function to get the description
  inline const xtool_ff_description& get_description(const xtool_ff a) {
    assert(XTOOL_FF_DICTIONARY[static_cast<int>(a)].value == a);
    return XTOOL_FF_DICTIONARY[static_cast<int>(a)];
  }

} // namespace mudock
