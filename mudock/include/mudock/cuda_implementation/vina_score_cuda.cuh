#pragma once

#include <mudock/compute/vina_score.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>

namespace mudock {
  template<>
  batch_multiple get_vina_score_batch_multiple<queue_cuda>(const int, std::shared_ptr<queue_cuda>);

  template<>
  void vina_score_kernel<queue_cuda>::operator()();
} // namespace mudock
