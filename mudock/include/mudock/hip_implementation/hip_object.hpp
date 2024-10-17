#pragma once

#include <cstdint>
#include <hip/hip_runtime.h>
#include <mudock/hip_implementation/hip_check_error_macro.hpp>

namespace mudock {

  template<class T>
  class hip_object {
    T* dev_ptr       = nullptr;
    std::size_t size = 0;
    const hipStream_t& stream;

  public:
    hip_object(const hipStream_t& _stream): stream(_stream) {}
    hip_object(const hip_object&) = delete;
    inline hip_object(hip_object&& other): stream(other.stream) {
      dev_ptr       = other.dev_ptr;
      size          = other.size;
      other.dev_ptr = nullptr;
      other.size    = 0;
    }
    ~hip_object() noexcept(false) {
      if (dev_ptr != nullptr)
        MUDOCK_CHECK(hipFreeAsync(dev_ptr, stream));
    }
    hip_object& operator=(const hip_object&) = delete;
    hip_object& operator=(hip_object&&)      = default;

    inline void alloc(const size_t num_elements) {
      if (size < num_elements) {
        if (dev_ptr != nullptr)
          MUDOCK_CHECK(hipFreeAsync(dev_ptr, stream));
        MUDOCK_CHECK(hipMallocAsync((void**)&dev_ptr, sizeof(T) * num_elements, stream));
      }
      size = num_elements;
    }
    // HIP set the const byte value https://rocm.docs.amd.com/projects/hipfort/en/docs-5.2.3/doxygen/html/interfacehipfort_1_1hipmemset.html
    inline void set_to_value(const int value) {
      MUDOCK_CHECK(hipMemsetAsync(dev_ptr, value, sizeof(T) * size, stream));
    }

    inline void copy_host2device(const T* const host) {
      MUDOCK_CHECK(hipMemcpyAsync(dev_ptr, host, sizeof(T) * size, hipMemcpyHostToDevice, stream));
    }
    inline void copy_device2host(T* const host) const {
      MUDOCK_CHECK(hipMemcpyAsync(host, dev_ptr, sizeof(T) * size, hipMemcpyDeviceToHost, stream));
    }

    [[nodiscard]] inline auto dev_pointer() const { return dev_ptr; }
    [[nodiscard]] inline auto num_elements() const { return size; }
    [[nodiscard]] inline auto get_stream() const { return stream; }
  };

} // namespace mudock
