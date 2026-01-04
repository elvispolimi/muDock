#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/genetic.hpp>
#include <mudock/hip_implementation/queue_hip.hpp>

namespace mudock {
  template<>
  void genetic_kernel<queue_hip>::operator()();
  template<>
  void genetic_kernel<queue_hip>::initialize();
  template<>
  void genetic_kernel<queue_hip>::finalize();
} // namespace mudock
