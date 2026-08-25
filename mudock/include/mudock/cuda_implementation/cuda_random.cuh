#pragma once

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <mudock/compute/object.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>

namespace mudock {
  struct cuda_random_object {
    cuda_random_object(std::shared_ptr<queue_cuda> q_): q(q_), state(q) {};
    cuda_random_object(const cuda_random_object &)            = delete;
    cuda_random_object(cuda_random_object &&)                 = delete;
    cuda_random_object &operator=(const cuda_random_object &) = delete;
    cuda_random_object &operator=(cuda_random_object &&)      = delete;

    void alloc(const std::size_t num_elements, const std::size_t seed);
    void alloc(const std::size_t num_elements);

    [[nodiscard]] inline auto dev_pointer() const { return state.dev_pointer(); }
    [[nodiscard]] inline curandStatePhilox4_32_10_t **dev_pointer_ref() { return state.dev_pointer_ref(); }

  private:
    std::shared_ptr<queue_cuda> q;
    object<curandStatePhilox4_32_10_t> state;
  };

} // namespace mudock
