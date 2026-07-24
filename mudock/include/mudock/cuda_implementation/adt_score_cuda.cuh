#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>

namespace mudock {
  template<>
  batch_multiple get_adt_score_batch_multiple<queue_cuda>(const int, std::shared_ptr<queue_cuda>);

  template<>
  void adt_score_kernel<queue_cuda>::operator()();
} // namespace mudock
