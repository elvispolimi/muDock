#pragma once

#include <limits>
#include <math.h>
#include <mudock/type_alias.hpp>
#include <array>

namespace mudock {
  static constexpr fp_type lambda{0.003627};
  static constexpr fp_type epsilon0{78.4};
  static constexpr fp_type A{-8.5525};
  static constexpr fp_type B = epsilon0 - A;
  static constexpr fp_type rk{7.7839};
  static constexpr fp_type lambda_B = -lambda * B;
  static constexpr fp_type min_epsilon{std::numeric_limits<fp_type>::epsilon()};

  inline fp_type calc_ddd_Mehler_Solmajer(const fp_type& distance) {
    fp_type epsilon = A + B / (fp_type{1} + rk * std::exp(lambda_B * distance));

    if (epsilon < min_epsilon) [[unlikely]] {
      epsilon = fp_type{1.0};
    }
    return epsilon;
  }

  static constexpr auto num_radius_tick     = std::size_t{2048};
  static constexpr auto num_radius_angstrom = fp_type{20.48};

  inline  auto compute_dielectric_ewds() {
    std::array<fp_type, num_radius_tick> result;
    result[0] = fp_type{1};
    for (std::size_t radius_index = 0; radius_index < num_radius_tick; ++radius_index) {
      const auto radius =
          static_cast<fp_type>(radius_index) * (num_radius_angstrom / static_cast<fp_type>(num_radius_tick));
      result[radius_index] =
          fp_type{1} / (radius * (A + B / (fp_type{1} + rk * std::exp(lambda_B * radius))));
      // result[radius_index] = fp_type{332} / (A + B / (fp_type{1} + rk * std::exp(lambda_B * radius)));
    }
    return result;
  }

} // namespace mudock
