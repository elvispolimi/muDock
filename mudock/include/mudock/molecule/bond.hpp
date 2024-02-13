#pragma once

#include <mudock/chem/bond_types.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  struct bond {
    index_type source;
    index_type dest;
    bond_type type;
  };

} // namespace mudock
