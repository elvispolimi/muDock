#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/adt_score_args.hpp>
#include <mudock/gh_implementation/queue_gh.hpp>

namespace mudock {
  void score_batch_gh(const adt_score_args& args);

  template<>
  void adt_score<queue_gh>::operator()();
} // namespace mudock
