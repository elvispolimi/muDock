#pragma once

#include <cstdint>
#include <mudock/chem/bond_types.hpp>

namespace mudock {

  struct bond {
    std::size_t source = std::size_t{0};
    std::size_t dest   = std::size_t{0};
    bond_type type     = bond_type::SINGLE;
    bool can_rotate    = false;
  };

} // namespace mudock
