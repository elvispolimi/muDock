#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/adt_score_kernel.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  template<>
  inline batch_multiple get_adt_score_batch_multiple<queue_cpp>(const int, std::shared_ptr<queue_cpp>) {
    return {10, 1};
  }

  template<>
  void adt_score_kernel<queue_cpp>::operator()();
} // namespace mudock
