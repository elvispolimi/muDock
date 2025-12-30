#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/adt_score_kernel.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  template<>
  inline int get_adt_score_batch<queue_cpp>(const int) {
    return 10;
  }

  template<>
  void adt_score_kernel<queue_cpp>::operator()();
} // namespace mudock
