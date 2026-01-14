#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/molecule.hpp>
#include <mudock/sycl_implementation/queue_sycl.hpp>

namespace mudock {
  template<>
  void geom_kernel<queue_sycl>::operator()();
} // namespace mudock
