#pragma once

#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/compute/geometric_transform.hpp>

namespace mudock {
  template<>
  void geom_kernel<queue_alpaka>::operator()();
} // namespace mudock
