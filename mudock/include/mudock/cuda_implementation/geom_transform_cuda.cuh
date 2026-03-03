#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>
#include <mudock/molecule.hpp>

namespace mudock {
  template<>
  batch_multiple get_geom_transform_batch_multiple<queue_cuda>(const int, std::shared_ptr<queue_cuda>);

  template<>
  void geom_kernel<queue_cuda>::operator()();
} // namespace mudock
