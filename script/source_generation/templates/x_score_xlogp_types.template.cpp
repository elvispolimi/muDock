#include <mudock/chem/x_score_xlogp_types.hpp>

namespace mudock {

  const std::array<xlogp_ff_description, {@ num_elements @}> XLOGP_FF_DICTIONARY = {{
{% for element in data %}
    {
      xlogp_ff::{@ element.name @},
      "{@ element.name @}",
      "{@ element.hbond @}",
      {@ element.hydrophobic_scale @}
    }{@ "," if not loop.last @}
{% endfor %}
  }};

} // namespace mudock
