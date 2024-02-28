#pragma once

#include <mudock/devices/context.hpp>
#include <pthread.h>
#include <thread>

namespace mudock {
  template<>
  device_context<language_types::CPP>::device_context(const index_type);

  template<>
  device_context<language_types::CPP>::~device_context();

} // namespace mudock
