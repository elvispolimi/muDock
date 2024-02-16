#pragma once

#include <mudock/devices/language_types.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  template<language_types l_t>
  class device_context {
  public:
    device_context(const index_type device_id);
    ~device_context();

    device_context(const device_context&)            = delete;
    device_context(device_context&&)                 = default;
    device_context& operator=(const device_context&) = delete;
    device_context& operator=(device_context&&)      = delete;
  };
} // namespace mudock
