#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/gh_implementation/queue_gh.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<>
  inline batch_multiple get_geom_transform_batch_multiple<queue_gh>(const int, std::shared_ptr<queue_gh>) {
    return {10, 1};
  }

  template<>
  void geom_kernel<queue_gh>::operator()();
} // namespace mudock
