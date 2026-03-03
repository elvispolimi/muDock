#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/molecule.hpp>
#include <mudock/sycl_implementation/queue_sycl.hpp>

namespace mudock {
  template<>
  batch_multiple get_geom_transform_batch_multiple<queue_sycl>(const int, std::shared_ptr<queue_sycl>);

  template<>
  void geom_kernel<queue_sycl>::operator()();
} // namespace mudock
