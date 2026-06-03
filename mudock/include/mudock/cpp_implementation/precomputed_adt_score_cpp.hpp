#pragma once

#include <mudock/compute/precomputed_adt_score.hpp>
#include <mudock/compute/precomputed_adt_score_kernel.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  template<>
  void precomputed_adt_score_kernel<queue_cpp>::operator()();
} // namespace mudock