#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/sycl_implementation/queue_sycl.hpp>

namespace mudock {
  template<>
  int get_adt_score_batch<queue_sycl>(const int, std::shared_ptr<queue_sycl>, const size_t);

  template<>
  void adt_score_kernel<queue_sycl>::operator()();
} // namespace mudock
