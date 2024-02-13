#include <algorithm>
#include <mudock/chem/elements.hpp>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/periodic_table.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  std::optional<element> parse(const std::string_view symbol) {
    const auto element_it = std::find_if(std::begin(ELEMENT_DICTIONARY),
                                         std::end(ELEMENT_DICTIONARY),
                                         [&symbol](const auto& e) { return e.symbol == symbol; });
    return element_it != std::end(ELEMENT_DICTIONARY) ? std::optional{element_it->value}
                                                      : std::optional<element>{};
  }

  const std::array<element_description, {@ num_elements @}> ELEMENT_DICTIONARY = {{{% for element in data %}
    {element::{@ element.value @}, "{@ element.symbol @}", "{@ element.name @}", {@ element.number @}, {@ element.valence @}},{% endfor %}
  }};

} // namespace mudock
