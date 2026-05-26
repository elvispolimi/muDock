#pragma once

#include <mudock/compute/adadelta.hpp>
#include <mudock/compute/adadelta_kernel.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  template<>
  void adadelta_kernel<queue_cpp>::operator()();
  
  template<>
  void adadelta_kernel<queue_cpp>::compute_gradients();

  template<>
  void adadelta_kernel<queue_cpp>::apply_adadelta();
} // namespace mudock
