#include <mudock/chem/autodock_types.hpp>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/autodock_types.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {
  const std::array<autodock_ff_description, {@ num_elements @}> AUTODOCK_FF_DICTIONARY = {{{%for element in data %}
    {autodock_ff::{@ element.value @}, "{@ element.name @}", {@ element.Rii |default(0.0) @}, {@ element.epsii |default(0.0) @}, {@ element.vol |default(0.0) @}, {@ element.solpar |default(0.0) @},  {@ element.Rij_hb |default(0.0) @},  {@ element.epsij_hb |default(0.0) @},  {@ element.hbond |default(0) @},},{% endfor %}
  }};

} // namespace mudock
