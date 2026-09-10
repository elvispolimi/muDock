#pragma once

#include <mudock/grid.hpp>
#include <mudock/type_alias.hpp>
#include <numeric>
#include <span>
#include <cmath>

namespace mudock {

  inline point<fp_type, 3> compute_center_of_mass(const std::span<const fp_type> x,
                                                  const std::span<const fp_type> y,
                                                  const std::span<const fp_type> z) {
    assert((x.size() == y.size()) && (y.size() == z.size()));
    return {std::accumulate(std::begin(x), std::end(x), fp_type{0}) / static_cast<fp_type>(x.size()),
            std::accumulate(std::begin(y), std::end(y), fp_type{0}) / static_cast<fp_type>(y.size()),
            std::accumulate(std::begin(z), std::end(z), fp_type{0}) / static_cast<fp_type>(z.size())};
  }

  inline point<fp_type, 3>
      compute_center_of_mass(const fp_type* x, const fp_type* y, const fp_type* z, const int num_atoms) {
    return {std::accumulate(x, x + num_atoms, fp_type{0}) / static_cast<fp_type>(num_atoms),
            std::accumulate(y, y + num_atoms, fp_type{0}) / static_cast<fp_type>(num_atoms),
            std::accumulate(z, z + num_atoms, fp_type{0}) / static_cast<fp_type>(num_atoms)};
  }

  inline fp_type
      compute_rmsd(const fp_type* x_A,
                   const fp_type* y_A,
                   const fp_type* z_A,
                   const fp_type* x_B,
                   const fp_type* y_B,
                   const fp_type* z_B,
                   const int num_atoms) {
    fp_type sum_sq_dist = 0.0;

    for (int i = 0; i < num_atoms; ++i) {
      const fp_type dx = x_A[i] - x_B[i];
      const fp_type dy = y_A[i] - y_B[i];
      const fp_type dz = z_A[i] - z_B[i];
      sum_sq_dist += dx * dx + dy * dy + dz * dz;
    }

    return std::sqrt(sum_sq_dist / static_cast<fp_type>(num_atoms));
  }

} // namespace mudock
