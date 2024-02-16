#pragma once

#include <mudock/devices/language_types.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  template<language_types k_t, typename T>
  class device_buffer {
  public:
    device_buffer(const index_type capacity): capacity(capacity){};
    ~device_buffer(){};

    device_buffer(const device_buffer&)            = delete;
    device_buffer(device_buffer&&)                 = delete;
    device_buffer& operator=(const device_buffer&) = delete;
    device_buffer& operator=(device_buffer&&)      = delete;

    virtual void copy_in(const T*, const index_type size);
    virtual void copy_out(T*);

  private:
    T* data_device;
    const index_type capacity;
    index_type size{0};
  };
} // namespace mudock
