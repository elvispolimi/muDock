#include <mudock/chem/bond_types.hpp>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/periodic_table.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  const std::array<bond_type_description, 5> BOND_DICTIONARY = {{
      {bond_type::SINGLE, "Single"},
      {bond_type::DOUBLE, "Double"},
      {bond_type::TRIPLE, "Triple"},
      {bond_type::AMIDE, "Amide"},
      {bond_type::AROMATIC, "Aromatic"},
  }};

} // namespace mudock
