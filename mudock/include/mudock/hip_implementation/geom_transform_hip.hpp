#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/hip_implementation/queue_hip.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<>
  batch_multiple get_geom_transform_batch_multiple<queue_hip>(const int, std::shared_ptr<queue_hip>);

  template<>
  void geom_kernel<queue_hip>::operator()();
} // namespace mudock
