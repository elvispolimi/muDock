#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/xsimd_implementation/queue_xsimd.hpp>

namespace mudock {

  template<>
  void adt_score_kernel<queue_xsimd>::operator()();
} // namespace mudock
