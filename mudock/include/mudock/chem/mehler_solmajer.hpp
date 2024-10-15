#pragma once

#include <limits>
#include <math.h>
#include <mudock/type_alias.hpp>

namespace mudock {
  fp_type calc_ddd_Mehler_Solmajer(fp_type distance) {
    const fp_type lambda{0.003627};
    const fp_type epsilon0{78.4};
    const fp_type A{-8.5525};
    const fp_type B = epsilon0 - A;
    const fp_type rk{7.7839};
    const fp_type lambda_B = -lambda * B;

    fp_type epsilon = A + B / (fp_type{1} + rk * expf(lambda_B * distance));

    if (epsilon < std::numeric_limits<fp_type>::epsilon()) {
      epsilon = fp_type{1.0};
    }
    return epsilon;
  }

} // namespace mudock
