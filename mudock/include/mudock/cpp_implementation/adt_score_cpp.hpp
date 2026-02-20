#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/bucket_size.hpp>
#include <mudock/compute/adt_score_kernel.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  template<>
  inline int get_adt_score_batch_multiple<queue_cpp>(const int, std::shared_ptr<queue_cpp>) {
    return 10;
  }

  template<>
  inline int get_adt_score_batch<queue_cpp>(const int atoms,
                                            std::shared_ptr<queue_cpp> q_b,
                                            const size_t max_bucket_size) {
    return resolve_bucket_size("CPP", atoms, max_bucket_size, [&]() {
      return get_adt_score_batch_multiple<queue_cpp>(atoms, q_b);
    });
  }

  template<>
  void adt_score_kernel<queue_cpp>::operator()();
} // namespace mudock
