#pragma once

#include <cstdint>
#include <cuda_runtime.h>

namespace mudock {

  template<class T>
  class cuda_object {
    T* dev_ptr       = nullptr;
    std::size_t size = 0;

  public:
    cuda_object()                   = default;
    cuda_object(const cuda_object&) = delete;
    inline cuda_object(cuda_object&& other) {
      dev_ptr       = other.dev_ptr;
      size          = other.size;
      other.dev_ptr = nullptr;
      other.size    = 0;
    }
    ~cuda_object() {
      if (dev_ptr != nullptr)
        CHECK(cudaFree(dev_ptr));
    }
    cuda cuda_object& operator=(const cuda_object&) = delete;
    cuda_object& operator=(cuda_object&&)           = default;

    inline void alloc(const size_t num_elements) {
      if (size < num_elements) {
        if (dev_ptr != nullptr)
          CHECK(cudaFree(&dev_ptr));
        CHECK(cudaMalloc(&dev_ptr, sizeof(T) * num_elements));
      }
      size = num_elements;
    }
    inline void set_to_value(const T value) { CHECK(cudaMemset(dev_ptr, value, size)); }

    inline void copy_host2device(const T* const host) {
      CHECK(cudaMemcpy(dev_ptr, host, sizeof(T) * size, cudaMemcpyHostToDevice))
    };
    inline void copy_odevice2host(T* const host) const {
      CHECK(cudaMemcpy(host, dev_ptr, sizeof(T) * size, cudaMemcpyDeviceToHost))
    };

    [[nodiscard]] inline auto dev_pointer() const { return dev_ptr; }
    [[nodiscard]] inline auto num_elements() const { return size; }
  };

} // namespace mudock
