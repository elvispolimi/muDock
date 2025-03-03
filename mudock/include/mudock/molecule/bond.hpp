#pragma once

#include <cstdint>
#include <mudock/chem/bond_types.hpp>

namespace mudock {

  struct bond {
    int source      = int{0};
    int dest        = int{0};
    bond_type type  = bond_type::SINGLE;
    bool can_rotate = false;
  };

} // namespace mudock
