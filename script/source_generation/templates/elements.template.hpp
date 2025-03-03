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

  // this is the list of all the known atoms
  enum class element : int {{% for element in data %}
    {@ element.value @} = {@ element.index @}, // {@ element.name @}{% endfor %}
  };

  // this is knowledge that we have about all the elements (that we need at least)
  struct element_description {
    element value;
    std::string_view symbol;
    std::string_view name;
    int number;
    int valence;
  };
  extern const std::array<element_description, {@num_elements@}> ELEMENT_DICTIONARY;

  // utility functions to work with them
  inline const element_description& get_description(const element e) {
    assert(ELEMENT_DICTIONARY[static_cast<int>(e)].value == e);
    return ELEMENT_DICTIONARY[static_cast<int>(e)];
  }
  std::optional<element> parse(const std::string_view symbol);

} // namespace mudock
