#pragma once

#include <cmath>
#include <mudock/type_alias.hpp>

namespace mudock {
  struct math {
    static constexpr auto pi        = fp_type{3.141592653589793238462643383279502884197};
    static constexpr auto pi_halved = pi / 2;
  };

  constexpr auto rad_to_deg(const fp_type rad) { return rad * fp_type{180} / math::pi; }
  constexpr auto deg_to_rad(const fp_type deg) { return deg * math::pi / fp_type{180}; }
} // namespace mudock
