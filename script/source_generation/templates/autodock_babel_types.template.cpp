#include <algorithm>
#include <mudock/chem/autodock_babel_types.hpp>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/autodock_babel_types.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {
  const std::array<autodock_babel_ff_description, {@ num_elements @}> AUTODOCK_BABEL_FF_DICTIONARY = {{{% for element in data %}
    {autodock_babel_ff::{@ element.value @}, "{@ element.name @}", element::{@ element.element @}},{% endfor %}
  }};

} // namespace mudock
