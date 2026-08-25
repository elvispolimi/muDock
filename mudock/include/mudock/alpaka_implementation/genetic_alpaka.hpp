#pragma once

#include <mudock/alpaka_implementation/buffer_alpaka.hpp>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/compute/genetic.hpp>

namespace mudock {
  template<>
  void genetic_kernel<queue_alpaka>::initialize();

  template<>
  void genetic_kernel<queue_alpaka>::operator()();

  template<>
  void genetic_kernel<queue_alpaka>::finalize();
} // namespace mudock
