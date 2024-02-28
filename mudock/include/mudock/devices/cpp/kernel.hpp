#pragma once

#include "mudock/devices/scratch.hpp"
#include <mudock/devices/context.hpp>
#include <mudock/devices/cpp/scratch.hpp>
#include <mudock/devices/kernel.hpp>

namespace mudock {
  template<>
  void dock_and_score<language_types::CPP>(language_scratchpad_impl<language_types::CPP>&,
                                           const device_context<language_types::CPP>&);
} // namespace mudock
