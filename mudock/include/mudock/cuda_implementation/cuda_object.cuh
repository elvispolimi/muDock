#pragma once

#include <cstdint>

namespace mudock {

  template<class T>
  class cuda_object {
    T* dev_ptr       = nullptr;
    std::size_t size = 0;

  public:
    cuda_object()                   = default;
    cuda_object(const cuda_object&) = delete;
    cuda_object(cuda_object&& other);
// TODO seems not supported by Polygeist
#ifdef MUDOCK_ENABLE_POLY
    ~cuda_object() = default;
#else
    ~cuda_object() noexcept(false);
#endif
    cuda_object& operator=(const cuda_object&) = delete;
    cuda_object& operator=(cuda_object&&)      = default;

    void alloc(const size_t num_elements);
    // CUDA memset set const byte value https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY.html#group__CUDART__MEMORY_1gf7338650f7683c51ee26aadc6973c63a
    void set_to_value(const int value);

    void copy_host2device(const T* const host);
    void copy_device2host(T* const host) const;

    [[nodiscard]] T* dev_pointer() const;
    [[nodiscard]] std::size_t num_elements() const;
  };
} // namespace mudock
