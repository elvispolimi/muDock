#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/genetic.hpp>
#include <mudock/sycl_implementation/queue_sycl.hpp>

namespace mudock {
  template<>
  void genetic_kernel<queue_sycl>::operator()();
  template<>
  void genetic_kernel<queue_sycl>::initialize();
  template<>
  void genetic_kernel<queue_sycl>::finalize();
} // namespace mudock
