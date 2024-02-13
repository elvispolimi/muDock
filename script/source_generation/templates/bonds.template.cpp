#include <mudock/chem/bond_types.hpp>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/periodic_table.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  const std::array<bond_type_description, {@num_bonds@}> BOND_DICTIONARY = {{{% for bond in data %}
    {bond_type::{@ bond.value @}, "{@ bond.name @}"},{% endfor %}
  }};

} // namespace mudock
