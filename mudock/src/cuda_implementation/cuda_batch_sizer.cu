#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cuda_implementation/cuda_batch_sizer.cuh>
#include <mudock/cuda_implementation/evaluate_fitness.cuh>
#include <mudock/utils.hpp>
#include <stdexcept>
#include <stdio.h>

#define BUCKET_MULTIPLIER 3

namespace mudock {
  int compute_batch_size(const int num_atoms, const int num_rotamers) {
    // TODO something useful
    return 1000;
  }
} // namespace mudock
