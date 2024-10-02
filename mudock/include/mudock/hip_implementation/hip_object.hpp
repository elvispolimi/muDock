#pragma once

#include <cstdint>
#include <hip/hip_runtime.h>
#include <mudock/hip_implementation/hip_check_error_macro.hpp>

namespace mudock {

  template<class T>
  class hip_object {
    T* dev_ptr       = nullptr;
    std::size_t size = 0;

  public:
    hip_object()                  = default;
    hip_object(const hip_object&) = delete;
    inline hip_object(hip_object&& other) {
      dev_ptr       = other.dev_ptr;
      size          = other.size;
      other.dev_ptr = nullptr;
      other.size    = 0;
    }
    ~hip_object() noexcept(false) {
      if (dev_ptr != nullptr)
        MUDOCK_CHECK(hipFree(dev_ptr));
    }
    hip_object& operator=(const hip_object&) = delete;
    hip_object& operator=(hip_object&&)      = default;

    inline void alloc(const size_t num_elements) {
      if (size < num_elements) {
        if (dev_ptr != nullptr)
          MUDOCK_CHECK(hipFree(dev_ptr));
        MUDOCK_CHECK(hipMalloc(&dev_ptr, sizeof(T) * num_elements));
      }
      size = num_elements;
    }
    inline void set_to_value(const T value) { MUDOCK_CHECK(hipMemset(dev_ptr, value, sizeof(T) * size)); }

    inline void copy_host2device(const T* const host) {
      MUDOCK_CHECK(hipMemcpy(dev_ptr, host, sizeof(T) * size, hipMemcpyHostToDevice));
    }
    inline void copy_device2host(T* const host) const {
      MUDOCK_CHECK(hipMemcpy(host, dev_ptr, sizeof(T) * size, hipMemcpyDeviceToHost));
    }

    [[nodiscard]] inline auto dev_pointer() const { return dev_ptr; }
    [[nodiscard]] inline auto num_elements() const { return size; }
  };

} // namespace mudock
