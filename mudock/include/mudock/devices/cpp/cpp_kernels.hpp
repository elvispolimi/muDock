#pragma once

#include <mudock/devices/device_contexts.hpp>
#include <mudock/devices/cpp/cpp_molecules.hpp>
#include <mudock/devices/device_kernels.hpp>

namespace mudock {
  template<>
  void dock_and_score<language_types::CPP>(molecules_scratchpad<language_types::CPP>&,
                                           const device_context<language_types::CPP>&);
} // namespace mudock
