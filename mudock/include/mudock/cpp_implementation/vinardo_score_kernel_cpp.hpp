#pragma once

#include <mudock/compute/vinardo_score_kernel.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  template<>
  void vinardo_score_kernel<queue_cpp>::operator()();
}
