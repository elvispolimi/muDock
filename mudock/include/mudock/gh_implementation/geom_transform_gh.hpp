#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/gh_implementation/queue_gh.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<>
  void geom_kernel<queue_gh>::operator()();
} // namespace mudock
