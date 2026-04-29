#pragma once

#include <limits>
#include <math.h>
#include <mudock/type_alias.hpp>

namespace mudock {
  namespace mehler_solmajer {
    static constexpr fp_type lambda      = static_cast<fp_type>(0.003627);
    static constexpr fp_type epsilon0    = static_cast<fp_type>(78.4);
    static constexpr fp_type A           = static_cast<fp_type>(-8.5525);
    static constexpr fp_type B           = epsilon0 - A;
    static constexpr fp_type rk          = static_cast<fp_type>(7.7839);
    static constexpr fp_type lambda_B    = -lambda * B;
    static constexpr fp_type min_epsilon = static_cast<fp_type>(std::numeric_limits<fp_type>::epsilon());
  } // namespace mehler_solmajer

  inline fp_type calc_ddd_Mehler_Solmajer(const fp_type& distance) {
    fp_type epsilon = mehler_solmajer::A +
                      mehler_solmajer::B /
                          (fp_type{1} + mehler_solmajer::rk * std::exp(mehler_solmajer::lambda_B * distance));

    // TODO
    //if (epsilon < mehler_solmajer::min_epsilon) [[unlikely]] {
    //  epsilon = fp_type{1.0};
    //}
    // return std::max(epsilon, fp_type{1.0});
    return epsilon;
  }
} // namespace mudock
