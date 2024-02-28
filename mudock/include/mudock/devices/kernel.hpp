#pragma once

#include <mudock/devices/context.hpp>
#include <mudock/devices/scratch.hpp>

namespace mudock {
  // This method has to implement the docking and the scoring functionalities
  // The method should be different based on the language type
  template<language_types l_t>
  void dock_and_score(language_scratchpad_impl<l_t>&, const device_context<l_t>&);
} // namespace mudock
