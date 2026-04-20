#pragma once

#include <mudock/compute/adt_quant_score.hpp>
#include <mudock/compute/adt_quant_score_kernel.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  template<>
  inline int get_adt_quant_score_batch<queue_cpp>(const int, std::shared_ptr<queue_cpp>) {
    return 10;
  }

  template<>
  void adt_quant_score_kernel<queue_cpp>::operator()();
} // namespace mudock
