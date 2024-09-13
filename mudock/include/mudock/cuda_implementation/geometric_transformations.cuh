#pragma once

#include <mudock/grid/point3D.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  __device__ void translate_molecule_cuda(fp_type* __restrict__ x,
                                          fp_type* __restrict__ y,
                                          fp_type* __restrict__ z,
                                          const fp_type* offset_x,
                                          const fp_type* offset_y,
                                          const fp_type* offset_z,
                                          const int num_atoms);

  __device__ void rotate_molecule_cuda(fp_type* __restrict__ x,
                                       fp_type* __restrict__ y,
                                       fp_type* __restrict__ z,
                                       const fp_type* angle_x,
                                       const fp_type* angle_y,
                                       const fp_type* angle_z,
                                       const int num_atoms);

  __device__ void rotate_fragment_cuda(fp_type* __restrict__ x,
                                       fp_type* __restrict__ y,
                                       fp_type* __restrict__ z,
                                       const int* bitmask,
                                       const int start_index,
                                       const int stop_index,
                                       const fp_type* angle,
                                       const int num_atoms);

} // namespace mudock
