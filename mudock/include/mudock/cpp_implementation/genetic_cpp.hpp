#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/genetic.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>
#include <random>

namespace mudock {
  template<>
  struct rand_state_type<queue_cpp> {
    using type = std::mt19937;
  };

  template<>
  void genetic_kernel<queue_cpp>::operator()();
  template<>
  void genetic_kernel<queue_cpp>::initialize();
  template<>
  void genetic_kernel<queue_cpp>::finalize();
} // namespace mudock
