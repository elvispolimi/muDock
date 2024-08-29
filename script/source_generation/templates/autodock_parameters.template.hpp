#pragma once

#include <mudock/type_alias.hpp>

namespace mudock {
  struct autodock_parameters {
    static constexpr fp_type coeff_vdW{{@ FE_coeff_vdW @}};
    static constexpr fp_type coeff_hbond{{@ FE_coeff_hbond @}};
    static constexpr fp_type coeff_estat{{@ FE_coeff_estat @}};
    static constexpr fp_type coeff_desolv{{@ FE_coeff_desolv @}};
    static constexpr fp_type coeff_tors{{@ FE_coeff_tors @}};
  };
} // namespace mudock
