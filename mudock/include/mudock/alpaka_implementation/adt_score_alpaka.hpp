#pragma once

#include <memory>
#include <mudock/alpaka_implementation/buffer_alpaka.hpp>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/compute/adt_score.hpp>

namespace mudock {
  template<>
  batch_multiple get_adt_score_batch_multiple<queue_alpaka>(const int, std::shared_ptr<queue_alpaka>);

  template<>
  void adt_score_kernel<queue_alpaka>::operator()();
} // namespace mudock
