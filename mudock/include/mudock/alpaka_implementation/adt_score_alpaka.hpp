#pragma once

#include <memory>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/compute/adt_score.hpp>

namespace mudock {
  template<>
  int get_adt_score_batch<queue_alpaka>(const int, std::shared_ptr<queue_alpaka>, const size_t);

  template<>
  void adt_score_kernel<queue_alpaka>::operator()();
} // namespace mudock
