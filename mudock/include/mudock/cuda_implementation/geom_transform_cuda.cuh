#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>
#include <mudock/molecule.hpp>

namespace mudock {
  template<>
  void geom_kernel<queue_cuda>::operator()();
} // namespace mudock
