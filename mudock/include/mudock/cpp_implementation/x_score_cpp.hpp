#pragma once

#include <mudock/compute/x_score.hpp>
#include <mudock/compute/x_score_kernel.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  template<>
  inline int get_x_score_batch<queue_cpp>(const int, std::shared_ptr<queue_cpp>, const size_t) {
    // put 1 stabilizes variance between MPI runs
    return 1;
  }

  template<>
  void x_score_kernel<queue_cpp>::operator()();
} // namespace mudock
