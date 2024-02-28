#pragma once

#include <cstddef>
#include <mudock/devices/kernel_types.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  template<kernel_type k_t, typename T>
  class device_buffer {
  public:
    device_buffer(const std::size_t capacity): capacity(capacity){};
    ~device_buffer(){};

    device_buffer(const device_buffer&)            = delete;
    device_buffer(device_buffer&&)                 = delete;
    device_buffer& operator=(const device_buffer&) = delete;
    device_buffer& operator=(device_buffer&&)      = delete;

    virtual void copy_in(const T* i_data, const std::size_t size);
    virtual void copy_out(T* o_data);

  private:
    T* data_device;
    const std::size_t capacity;
    std::size_t size{0};
  };
} // namespace mudock
