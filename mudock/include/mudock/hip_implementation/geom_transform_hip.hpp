#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/hip_implementation/queue_hip.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<>
  void geometric<queue_hip>::operator()();
} // namespace mudock
