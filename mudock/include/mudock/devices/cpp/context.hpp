#pragma once

#include <cstdint>
#include <mudock/devices/context.hpp>
#include <pthread.h>
#include <thread>

namespace mudock {
  template<>
  device_context<kernel_type::CPP>::device_context(const std::size_t);

  template<>
  device_context<kernel_type::CPP>::~device_context();

} // namespace mudock
