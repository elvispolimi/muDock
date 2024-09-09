#pragma once

#include <mudock/grid/point3D.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  __device__ void translate_molecule(fp_type* __restrict__ x,
                                     fp_type* __restrict__ y,
                                     fp_type* __restrict__ z,
                                     const fp_type* offset_x,
                                     const fp_type* offset_y,
                                     const fp_type* offset_z,
                                     const std::size_t num_atoms);

  __device__ void rotate_molecule(fp_type* __restrict__ x,
                                  fp_type* __restrict__ y,
                                  fp_type* __restrict__ z,
                                  const fp_type* angle_x,
                                  const fp_type* angle_y,
                                  const fp_type* angle_z,
                                  const std::size_t num_atoms);

  __device__ void rotate_fragment(fp_type* __restrict__ x,
                                  fp_type* __restrict__ y,
                                  fp_type* __restrict__ z,
                                  const typename fragments<static_containers>::value_type* bitmask,
                                  const std::size_t start_index,
                                  const std::size_t stop_index,
                                  const fp_type* angle,
                                  const std::size_t num_atoms);

} // namespace mudock
