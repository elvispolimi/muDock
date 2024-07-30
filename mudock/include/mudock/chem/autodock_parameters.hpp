#pragma once

#include <mudock/type_alias.hpp>

namespace mudock {
  struct autodock_parameters {
    static constexpr fp_type coeff_vdW{0.1560};
    static constexpr fp_type coeff_hbond{0.0974};
    static constexpr fp_type coeff_estat{0.1465};
    static constexpr fp_type coeff_desolv{0.1159};
    static constexpr fp_type coeff_tors{0.2744};
  };
} // namespace mudock
