#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/gh_implementation/queue_gh.hpp>

namespace mudock {

  template<>
  void adt_score_kernel<queue_gh>::operator()();
} // namespace mudock
