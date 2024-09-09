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
             const std::size_t* __restrict__ fragments_start_index,
             const std::size_t* __restrict__ fragments_stop_index,
             const std::size_t num_rotamers,
             const std::size_t stride_atoms,
             const std::size_t num_atoms);

} // namespace mudock
