#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/xsimd_implementation/queue_xsimd.hpp>

namespace mudock {
  template<>
  inline int get_adt_score_batch<queue_xsimd>(const int) {
    return 10;
  }
  template<>
  void adt_score_kernel<queue_xsimd>::operator()();
} // namespace mudock
