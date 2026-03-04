#pragma once

#include <mudock/compute/vina_score.hpp>
#include <mudock/compute/vina_score_kernel.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  template<>
  inline int get_vina_score_batch<queue_cpp>(const int, std::shared_ptr<queue_cpp>) {
    return 10;
  }

  template<>
  void vina_score_kernel<queue_cpp>::operator()();
}
