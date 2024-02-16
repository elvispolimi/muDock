#pragma once

#include <mudock/devices/device_contexts.hpp>
#include <mudock/devices/molecules.hpp>

namespace mudock {
  template<language_types l_t>
  void dock_and_score(molecules_scratchpad<l_t>&, const device_context<l_t>&);
} // namespace mudock
