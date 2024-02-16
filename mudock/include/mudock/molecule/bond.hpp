#pragma once

#include <mudock/chem/bond_types.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  struct bond {
    index_type source = index_type{0};
    index_type dest   = index_type{0};
    bond_type type    = bond_type::SINGLE;
    bool can_rotate   = false;
  };

} // namespace mudock
