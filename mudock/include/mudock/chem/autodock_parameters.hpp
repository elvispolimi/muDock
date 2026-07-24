#pragma once

#include <mudock/type_alias.hpp>

namespace mudock {
  struct autodock_parameters {
    // TODO hardcode from param_string_4_1
    static constexpr fp_type coeff_vdW    = static_cast<fp_type>(0.1662);
    static constexpr fp_type coeff_hbond  = static_cast<fp_type>(0.1209);
    static constexpr fp_type coeff_estat  = static_cast<fp_type>(0.1406);
    static constexpr fp_type coeff_desolv = static_cast<fp_type>(0.1322);
    static constexpr fp_type coeff_tors   = static_cast<fp_type>(0.2983);
  };
} // namespace mudock
