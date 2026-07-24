#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/hip_implementation/queue_hip.hpp>

namespace mudock {
  template<>
  batch_multiple get_adt_score_batch_multiple<queue_hip>(const int, std::shared_ptr<queue_hip>);

  template<>
  void adt_score_kernel<queue_hip>::operator()();
} // namespace mudock
