#pragma once

#include <limits>
#include <math.h>
#include <mudock/type_alias.hpp>
#include <array>

namespace mudock {
  namespace mehler_solmajer{
  static constexpr fp_type lambda{0.003627};
  static constexpr fp_type epsilon0{78.4};
  static constexpr fp_type A{-8.5525};
  static constexpr fp_type B = epsilon0 - A;
  static constexpr fp_type rk{7.7839};
  static constexpr fp_type lambda_B = -lambda * B;
  static constexpr fp_type min_epsilon{std::numeric_limits<fp_type>::epsilon()};
  }

  inline fp_type calc_ddd_Mehler_Solmajer(const fp_type& distance) {
    fp_type epsilon = mehler_solmajer::A + mehler_solmajer::B / (fp_type{1} + mehler_solmajer::rk * std::exp(mehler_solmajer::lambda_B * distance));

    if (epsilon < mehler_solmajer::min_epsilon) [[unlikely]] {
      epsilon = fp_type{1.0};
    }
    return epsilon;
  }
} // namespace mudock
