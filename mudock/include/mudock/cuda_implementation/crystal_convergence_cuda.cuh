#pragma once

#include <mudock/compute/crystal_convergence.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>

namespace mudock {
  template<>
  batch_multiple get_crystal_conv_batch_multiple<queue_cuda>(const int, std::shared_ptr<queue_cuda>);

  template<>
  void crystal_convergence_kernel<queue_cuda>::operator()();
} // namespace mudock
