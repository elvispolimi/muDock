#pragma once

#include "chromosome.hpp"

namespace mudock {
  void apply(fp_type* __restrict__ x,
             fp_type* __restrict__ y,
             fp_type* __restrict__ z,
             const chromosome& c,
             const int num_atoms,
             const int num_rotamers,
             const int* __restrict__ frag_masks,
             const int* __restrict__ frag_start_indexes,
             const int* __restrict__ frag_stop_indexes);
}
