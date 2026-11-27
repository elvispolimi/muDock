#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/molecule.hpp>
#include <mudock/xsimd_implementation/queue_xsimd.hpp>

namespace mudock {
  template<>
  void geometric<queue_xsimd>::operator()();
} // namespace mudock
