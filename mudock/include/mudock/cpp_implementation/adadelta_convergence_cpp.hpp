#pragma once

#include <mudock/compute/adadelta_convergence.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  template<>
  void adadelta_convergence_kernel<queue_cpp>::operator()();
} // namespace mudock
