#pragma once

#include <array>
#include <cuda.h>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <random>
#include <span>

namespace mudock {
  void apply_omp(fp_type* __restrict__ x,
                            fp_type* __restrict__ y,
                            fp_type* __restrict__ z,
                            const chromosome& chromosome,
                            const int* __restrict__ fragments,
                            const int* __restrict__ fragments_start_index,
                            const int* __restrict__ fragments_stop_index,
                            const int num_rotamers,
                            const int stride_atoms,
                            const int num_atoms);

} // namespace mudock
