#pragma once

#include <array>
#include <cuda.h>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cuda_implementation/geometric_transformations.cuh>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <random>
#include <span>

namespace mudock {

  __device__ inline void apply_cuda(fp_type* __restrict__ x,
                             fp_type* __restrict__ y,
                             fp_type* __restrict__ z,
                             const chromosome& chromosome,
                             const int* __restrict__ fragments,
                             const int* __restrict__ fragments_start_index,
                             const int* __restrict__ fragments_stop_index,
                             const int num_rotamers,
                             const int stride_atoms,
                             const int num_atoms) {
    // apply rigid transformations
    translate_molecule_cuda(x, y, z, &chromosome[0], &chromosome[1], &chromosome[2], num_atoms);
    rotate_molecule_cuda(x, y, z, &chromosome[3], &chromosome[4], &chromosome[5], num_atoms);

    // change the molecule shape
    for (int i = 0; i < num_rotamers; ++i) {
      const int* bitmask = fragments + i * stride_atoms;
      rotate_fragment_cuda(x,
                           y,
                           z,
                           bitmask,
                           fragments_start_index[i],
                           fragments_stop_index[i],
                           &chromosome[6 + i],
                           num_atoms);
    }
  }

} // namespace mudock
