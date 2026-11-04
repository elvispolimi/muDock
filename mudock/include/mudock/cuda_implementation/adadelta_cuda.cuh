#pragma once

#include <mudock/compute/adadelta.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>

namespace mudock {
  template<>
  batch_multiple get_adadelta_batch_multiple<queue_cuda>(const int, std::shared_ptr<queue_cuda>);

  template<>
  void adadelta_kernel<queue_cuda>::operator()();
} // namespace mudock
