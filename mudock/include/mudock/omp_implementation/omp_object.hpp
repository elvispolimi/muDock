#pragma once

#include <cstdint>
#include <omp.h>
#include <stdexcept>

namespace mudock {

  template<class T>
  class omp_object {
    T* dev_ptr       = nullptr;
    std::size_t size = 0;

  public:
    omp_object()                  = default;
    omp_object(const omp_object&) = delete;
    inline omp_object(omp_object&& other) {
      dev_ptr       = other.dev_ptr;
      size          = other.size;
      other.dev_ptr = nullptr;
      other.size    = 0;
    }
    ~omp_object() noexcept(false) {
      if (dev_ptr != nullptr)
        omp_target_free(dev_ptr, omp_get_default_device());
    }
    omp_object& operator=(const omp_object&) = delete;
    omp_object& operator=(omp_object&&)      = default;

    inline void alloc(const size_t num_elements) {
      if (size < num_elements) {
        if (dev_ptr != nullptr)
          omp_target_free(dev_ptr, omp_get_default_device());
        dev_ptr = (T*) omp_target_alloc(num_elements * sizeof(T), omp_get_default_device());
      }
      size = num_elements;
    }

    inline void set_to_value([[maybe_unused]] const int value) {
      static_assert(true, "OpenMP memset not defined");
    }

    inline void copy_host2device(const T* const host) {
      // Copy memory from host to device
      omp_target_memcpy(dev_ptr,
                        host,
                        size * sizeof(T),
                        0,
                        0,
                        omp_get_default_device(),
                        omp_get_initial_device());
    }
    inline void copy_device2host(T* const host) const {
      // Copy memory from host to device
      omp_target_memcpy(host,
                        dev_ptr,
                        size * sizeof(T),
                        0,
                        0,
                        omp_get_initial_device(),
                        omp_get_default_device());
    }

    [[nodiscard]] inline auto dev_pointer() const { return dev_ptr; }
    [[nodiscard]] inline auto num_elements() const { return size; }
  };

} // namespace mudock
