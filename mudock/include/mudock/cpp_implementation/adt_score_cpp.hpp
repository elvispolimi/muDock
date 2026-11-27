#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/adt_score_args.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  void score_batch_cpp(const adt_score_args& args);

  template<>
  void adt_score<queue_cpp>::operator()();
} // namespace mudock
