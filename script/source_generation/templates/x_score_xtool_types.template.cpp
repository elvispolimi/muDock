#include <mudock/chem/x_score_xtool_types.hpp>

namespace mudock {

  const std::array<xtool_ff_description, {@ num_elements @}> XTOOL_FF_DICTIONARY = {{
{% for element in data %}
    {
      xtool_ff::{@ element.name @},
      "{@ element.name @}",
      {@ element.atomic_weight @},
      {@ element.vdw_radius @},
      {@ element.vdw_potential @},
      {@ element.par_charge @},
      "{@ element.hbond @}"
    }{@ "," if not loop.last @}
{% endfor %}
  }};

} // namespace mudock
