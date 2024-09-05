#pragma once

#include <mudock/grid.hpp>
#include <mudock/type_alias.hpp>
#include <numeric>
#include <span>

namespace mudock {

  inline point3D compute_center_of_mass(const std::span<const fp_type> x,
                                        const std::span<const fp_type> y,
                                        const std::span<const fp_type> z) {
    return point3D{std::accumulate(std::begin(x), std::end(x), fp_type{0}) / static_cast<fp_type>(x.size()),
                   std::accumulate(std::begin(y), std::end(y), fp_type{0}) / static_cast<fp_type>(y.size()),
                   std::accumulate(std::begin(z), std::end(z), fp_type{0}) / static_cast<fp_type>(z.size())};
  }

} // namespace mudock
