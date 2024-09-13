#pragma once

#include <array>
#include <mudock/molecule.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/type_alias.hpp>
#include <cuda.h>
#include <random>
#include <span>

namespace mudock {
  __device__ void apply_cuda(fp_type* __restrict__ x,
             fp_type* __restrict__ y,
             fp_type* __restrict__ z,
             const fp_type* chromosome,
             const int* __restrict__ fragments,
             const int* __restrict__ fragments_start_index,
             const int* __restrict__ fragments_stop_index,
             const int num_rotamers,
             const int stride_atoms,
             const int num_atoms);

} // namespace mudock
