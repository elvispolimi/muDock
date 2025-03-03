#pragma once

#include <mudock/grid/point3D.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  void translate_molecule(std::span<fp_type> x,
                          std::span<fp_type> y,
                          std::span<fp_type> z,
                          const fp_type offset_x,
                          const fp_type offset_y,
                          const fp_type offset_z);

  void rotate_molecule(std::span<fp_type> x,
                       std::span<fp_type> y,
                       std::span<fp_type> z,
                       const fp_type angle_x,
                       const fp_type angle_y,
                       const fp_type angle_z);

  void rotate_fragment(std::span<fp_type> x,
                       std::span<fp_type> y,
                       std::span<fp_type> z,
                       std::span<const typename fragments<static_containers>::value_type> bitmask,
                       const int start_index,
                       const int stop_index,
                       const fp_type angle);

} // namespace mudock
