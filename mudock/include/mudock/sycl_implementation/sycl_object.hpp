#pragma once

#include <cstdint>
#include <sycl/sycl.hpp>

namespace mudock {

  template<class T>
  class sycl_object {
    sycl::queue& queue;
    T* dev_ptr       = nullptr;
    std::size_t size = 0;

  public:
    sycl_object(sycl::queue& _queue): queue(_queue){};
    sycl_object(const sycl_object&) = delete;
    inline sycl_object(sycl_object&& other): queue(other.queue) {
      // TODO the other remains with the same queue
      dev_ptr       = other.dev_ptr;
      size          = other.size;
      other.dev_ptr = nullptr;
      other.size    = 0;
    }
    ~sycl_object() noexcept(false) {
      if (dev_ptr != nullptr)
        sycl::free(dev_ptr, queue);
    }
    sycl_object& operator=(const sycl_object&) = delete;
    sycl_object& operator=(sycl_object&&)      = default;

    inline void alloc(const size_t num_elements) {
      if (size < num_elements) {
        if (dev_ptr != nullptr)
          sycl::free(dev_ptr, queue);
        dev_ptr = (T*) sycl::malloc_device(sizeof(T) * num_elements, queue);
      }
      size = num_elements;
    }
    // Set const byte value, value is interpreted as an unsigned char
    // https://registry.khronos.org/SYCL/specs/sycl-2020/html/sycl-2020.html
    // 4.9.4.3. SYCL functions for explicit memory operations
    inline void set_to_value(const int value) { queue.memset(dev_ptr, value, sizeof(T) * size).wait(); }

    inline void copy_host2device(const T* const host) {
      queue.memcpy(dev_ptr, host, sizeof(T) * size).wait();
    }
    inline void copy_device2host(T* const host) const {
      queue.memcpy(host, dev_ptr, sizeof(T) * size).wait();
    }

    [[nodiscard]] inline auto dev_pointer() const { return dev_ptr; }
    [[nodiscard]] inline auto num_elements() const { return size; }
  };

} // namespace mudock
