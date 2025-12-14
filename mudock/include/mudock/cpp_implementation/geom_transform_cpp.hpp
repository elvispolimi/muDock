#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<>
  void geom_kernel<queue_cpp>::operator()();
} // namespace mudock
