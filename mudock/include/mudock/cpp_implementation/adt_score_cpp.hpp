#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/adt_score_kernel.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  template<>
  void adt_score_kernel<queue_cpp>::operator()();
} // namespace mudock
