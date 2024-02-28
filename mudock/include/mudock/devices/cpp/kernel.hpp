#pragma once

#include <mudock/devices/context.hpp>
#include <mudock/devices/cpp/scratch.hpp>
#include <mudock/devices/kernel.hpp>

namespace mudock {
  template<>
  void dock_and_score<kernel_type::CPP>(kernel_scratchpad_impl<kernel_type::CPP>& d_scratchpad,
                                        const device_context<kernel_type::CPP>& device_c);
} // namespace mudock
