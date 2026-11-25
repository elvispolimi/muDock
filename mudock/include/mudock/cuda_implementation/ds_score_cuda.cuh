#pragma once

#include <mudock/compute/drugscore/ds_score.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>

namespace mudock {
  template<>
  void ds_score_kernel<queue_cuda>::operator()();
} // namespace mudock
