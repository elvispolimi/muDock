#pragma once

#include <mudock/compute/adadelta_kernel.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>

namespace mudock {
  template<>
  void adadelta_kernel<queue_cuda>::operator()();
} // namespace mudock
