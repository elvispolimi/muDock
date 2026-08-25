#pragma once

#include <array>
#include <cassert>
#include <mudock/type_alias.hpp>
#include <string_view>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/x_score_xlogp_types.json and subsequently modified manually
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  // List of all known atoms for X-LOGP force field
  enum class xlogp_ff : int {
{% for element in data %}
    {@ element.name @} = {@ element.index @}, // {@ element.name @}
{% endfor %}
  };

  // Knowledge about the X-LOGP force field parameters
  struct xlogp_ff_description {
    xlogp_ff value;
    std::string_view name;
    std::string_view hbond;
    fp_type hydrophobic_scale;
  };

  extern const std::array<xlogp_ff_description, {@ num_elements @}> XLOGP_FF_DICTIONARY;

  // Utility function to get the description
  inline const xlogp_ff_description& get_description(const xlogp_ff a) {
    assert(XLOGP_FF_DICTIONARY[static_cast<int>(a)].value == a);
    return XLOGP_FF_DICTIONARY[static_cast<int>(a)];
  }

} // namespace mudock
