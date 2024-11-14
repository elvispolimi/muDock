#pragma once

#include <limits>
#include <math.h>
#include <mudock/type_alias.hpp>
#include <array>

namespace mudock {
  static constexpr fp_type ms_lambda{0.003627};
  static constexpr fp_type ms_epsilon0{78.4};
  static constexpr fp_type ms_A{-8.5525};
  static constexpr fp_type ms_B = ms_epsilon0 - ms_A;
  static constexpr fp_type ms_rk{7.7839};
  static constexpr fp_type ms_lambda_B = -ms_lambda * ms_B;
  static constexpr fp_type ms_min_epsilon{std::numeric_limits<fp_type>::epsilon()};

  inline fp_type calc_ddd_Mehler_Solmajer(const fp_type& distance) {
    fp_type epsilon = ms_A + ms_B / (fp_type{1} + ms_rk * std::exp(ms_lambda_B * distance));

    if (epsilon < ms_min_epsilon) [[unlikely]] {
      epsilon = fp_type{1.0};
    }
    return epsilon;
  }
} // namespace mudock
