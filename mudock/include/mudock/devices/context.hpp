#pragma once

#include <cstdint>
#include <mudock/devices/kernel_types.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  // Each kernel requires a different way to handle the device's context
  // For ech new kernel type a new device context has to be provided
  // Each thread will have each own context, of the selected device
  template<kernel_type l_t>
  class device_context {
  public:
    // Each device context requires the ID of the target device
    device_context(const std::size_t device_id);
    ~device_context();

    // The device context has to be unique per (thread,device) pairs
    device_context(const device_context&)            = delete;
    device_context(device_context&&)                 = default;
    device_context& operator=(const device_context&) = delete;
    device_context& operator=(device_context&&)      = delete;
  };
} // namespace mudock
