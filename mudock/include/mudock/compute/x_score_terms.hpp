#pragma once

#include <mudock/type_alias.hpp>

namespace mudock {

  // XScore's calculated outputs: used to access the output of the kernel
  enum x_score_term : int {
    x_term_vdw = 0, // van der Waals
    x_term_hb  = 1, // hydrogen bond
    x_term_hp  = 2, // hydrophobic pairwise
    x_term_rt  = 3, // rotor penalty
    x_term_pkd = 4, // predicted -log(Kd) (regression over the four terms above)
    x_term_count = 5
  };

  // Regression coefficients of implemented XScore's terms: (vdw, hb, hp, rt; hm, hs not implemented)
  //
  //   pKd = c0 + cvdw*vdw + chb*hb + chp*hp + crt*rt          (XScore's pkd1)
  //
  static constexpr fp_type x_hpscore_cvdw = fp_type{0.004f};
  static constexpr fp_type x_hpscore_chb  = fp_type{0.054f};
  static constexpr fp_type x_hpscore_chp  = fp_type{0.009f};
  static constexpr fp_type x_hpscore_crt  = fp_type{-0.061f};
  static constexpr fp_type x_hpscore_c0   = fp_type{3.441f};

  // Regression function
  [[nodiscard]] constexpr fp_type compute_x_score_pkd(const fp_type vdw,
                                                      const fp_type hb,
                                                      const fp_type hp,
                                                      const fp_type rt) {
    return x_hpscore_c0 + x_hpscore_cvdw * vdw + x_hpscore_chb * hb + x_hpscore_chp * hp +
           x_hpscore_crt * rt;
  }

} // namespace mudock
