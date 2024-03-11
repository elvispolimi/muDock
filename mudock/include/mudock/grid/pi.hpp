#pragma once

#include <cmath>
#include <mudock/type_alias.hpp>

namespace mudock {
  struct math {
    static constexpr auto pi = coordinate_type{3.141592653589793238462643383279502884197};
  };

  constexpr auto rad_to_deg(const coordinate_type angle) { return angle * coordinate_type{180} / math::pi; }
} // namespace mudock
