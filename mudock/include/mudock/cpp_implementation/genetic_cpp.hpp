#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/genetic.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {

  template<>
  void genetic_kernel<queue_cpp>::operator()();
  template<>
  void genetic_kernel<queue_cpp>::initialize();
  template<>
  void genetic_kernel<queue_cpp>::finalize();
} // namespace mudock
