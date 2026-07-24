#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/xsimd_implementation/queue_xsimd.hpp>

namespace mudock {
  template<>
  inline batch_multiple get_adt_score_batch_multiple<queue_xsimd>(const int, std::shared_ptr<queue_xsimd>) {
    return {10, 1};
  }
  template<>
  void adt_score_kernel<queue_xsimd>::operator()();
} // namespace mudock
