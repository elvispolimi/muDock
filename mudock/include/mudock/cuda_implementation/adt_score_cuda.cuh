#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>

namespace mudock {
  template<>
  int get_adt_score_batch<queue_cuda>(const int);

  template<>
  void adt_score_kernel<queue_cuda>::operator()();
} // namespace mudock
