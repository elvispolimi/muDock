#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/genetic.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>

namespace mudock {
  template<>
  void genetic_kernel<queue_cuda>::operator()();
  template<>
  void genetic_kernel<queue_cuda>::initialize();
  template<>
  void genetic_kernel<queue_cuda>::finalize();
} // namespace mudock
