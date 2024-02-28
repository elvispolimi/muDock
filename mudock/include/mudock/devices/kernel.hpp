#pragma once

#include <mudock/devices/context.hpp>
#include <mudock/devices/scratch.hpp>

namespace mudock {
  // This method has to implement the docking and the scoring functionalities
  // The method should be different based on the language type
  template<kernel_type l_t>
  void dock_and_score(kernel_scratchpad_impl<l_t>& d_scratchpad, const device_context<l_t>& device_c);
} // namespace mudock
