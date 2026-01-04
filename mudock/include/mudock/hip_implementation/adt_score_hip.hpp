#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/hip_implementation/queue_hip.hpp>

namespace mudock {
  template<>
  int get_adt_score_batch<queue_hip>(const int);

  template<>
  void adt_score_kernel<queue_hip>::operator()();
} // namespace mudock
