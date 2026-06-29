#pragma once

#include <mudock/compute/adadelta.hpp>
#include <mudock/compute/adadelta_kernel.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>

namespace mudock {
  template<>
  void adadelta_kernel<queue_cuda>::operator()();
  
  template<>
  void adadelta_kernel<queue_cuda>::compute_gradients();

  template<>
  void adadelta_kernel<queue_cuda>::apply_adadelta();
} // namespace mudock
