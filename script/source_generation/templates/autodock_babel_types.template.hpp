#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <mudock/chem/elements.hpp>
#include <optional>
#include <string_view>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/autodock_babel_types.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  // this is the list of all the known atoms
  enum class autodock_babel_ff : int {{% for element in data %}
    {@ element.value @} = {@ element.index @}, // {@ element.value @}{% endfor %}
  };

  // this is knowledge tha/home/gianmarco/repo/muDock/mudock/src/compute/parse_ids.cpp:22:38:p/home/gianmarco/repo/muDock/mudock/src/compute/parse_ids.cpp:22:38:t we have about all the elements (that we need at least)
  struct autodock_babel_ff_description {
    autodock_babel_ff value;
    std::string_view name;
    element base_element;
  };
  extern const std::array < autodock_babel_ff_description, { @num_elements @ }
  > AUTODOCK_BABEL_FF_DICTIONARY;

  // utility functions to work with them
  inline const autodock_babel_ff_description& get_description(const autodock_babel_ff a) {
    assert(AUTODOCK_BABEL_FF_DICTIONARY[static_cast<int>(a)].value == a);
    return AUTODOCK_BABEL_FF_DICTIONARY[static_cast<int>(a)];
  }

} // namespace mudock
