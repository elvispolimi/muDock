#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/gh_implementation/queue_gh.hpp>

namespace mudock {
  template<>
  inline batch_multiple get_adt_score_batch_multiple<queue_gh>(const int, std::shared_ptr<queue_gh>) {
    return {10, 1};
  }
  template<>
  void adt_score_kernel<queue_gh>::operator()();
} // namespace mudock
