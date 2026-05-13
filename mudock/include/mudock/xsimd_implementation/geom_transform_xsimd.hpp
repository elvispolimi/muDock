#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/molecule.hpp>
#include <mudock/xsimd_implementation/queue_xsimd.hpp>

namespace mudock {
  template<>
  inline batch_multiple get_geom_transform_batch_multiple<queue_xsimd>(const int, std::shared_ptr<queue_xsimd>) {
    return {10, 1};
  }

  template<>
  void geom_kernel<queue_xsimd>::operator()();
} // namespace mudock
