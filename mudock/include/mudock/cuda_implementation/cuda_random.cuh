#pragma once

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <mudock/cuda_implementation/cuda_object.cuh>

namespace mudock {
  struct cuda_random_object: private cuda_object<curandState> {
    cuda_random_object()                                      = default;
    cuda_random_object(const cuda_random_object &)            = delete;
    cuda_random_object(cuda_random_object &&)                 = default;
    cuda_random_object &operator=(const cuda_random_object &) = delete;
    cuda_random_object &operator=(cuda_random_object &&)      = default;

    void alloc(const std::size_t num_elements);
    [[nodiscard]] inline auto dev_pointer() const { return cuda_object<curandState>::dev_pointer(); }
  };

} // namespace mudock