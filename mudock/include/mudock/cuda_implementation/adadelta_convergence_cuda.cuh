#pragma once

#include <mudock/compute/adadelta_convergence.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>

namespace mudock {
  template<>
  void adadelta_convergence_kernel<queue_cuda>::operator()(std::size_t completed_iterations,
                                                          bool finalize);
} // namespace mudock
