#include <mudock/chem/autodock_types.hpp>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/autodock_types.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {
const std::array<autodock_ff_description, {@ num_elements @}> AUTODOCK_FF_DICTIONARY = {{ {% for element in data %}
  {
    autodock_ff::{@ element.value @},
    "{@ element.name @}",
    static_cast<fp_type>({@ element.Rii | default(0.0) @}),
    static_cast<fp_type>({@ element.epsii | default(0.0) @}),
    static_cast<fp_type>({@ element.vol | default(0.0) @}),
    static_cast<fp_type>({@ element.solpar | default(0.0) @}),
    static_cast<fp_type>({@ element.Rij_hb | default(0.0) @}),
    static_cast<fp_type>({@ element.epsij_hb | default(0.0) @}),
    {@ element.hbond | default(0) @}
  }{% if not loop.last %},{% endif %}
{% endfor %} }};
} // namespace mudock
