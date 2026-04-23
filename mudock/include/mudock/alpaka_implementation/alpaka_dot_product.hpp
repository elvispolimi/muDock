#pragma once

#include <span>

namespace mudock::alpaka_demo {

  double dot_product(std::span<const float> lhs, std::span<const float> rhs);
  double dot_product_reduction(std::span<const float> lhs, std::span<const float> rhs);

} // namespace mudock::alpaka_demo
