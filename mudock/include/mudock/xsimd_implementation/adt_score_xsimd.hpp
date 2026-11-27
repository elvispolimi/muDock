#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/adt_score_args.hpp>
#include <mudock/xsimd_implementation/queue_xsimd.hpp>

namespace mudock {
  void score_batch_xsimd(const adt_score_args& args);

  template<>
  void adt_score<queue_xsimd>::operator()();
} // namespace mudock
