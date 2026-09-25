#pragma once

#include <mudock/compute/crystal_convergence.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  template<>
  void crystal_convergence_kernel<queue_cpp>::operator()();
} // namespace mudock
