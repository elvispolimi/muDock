#pragma once

#include <mudock/compute/adt_score.hpp>
#include <mudock/sycl_implementation/queue_sycl.hpp>

namespace mudock {
  template<>
  batch_multiple get_adt_score_batch_multiple<queue_sycl>(const int, std::shared_ptr<queue_sycl>);

  template<>
  void adt_score_kernel<queue_sycl>::operator()();
} // namespace mudock
