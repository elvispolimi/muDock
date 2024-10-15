#include <mudock/cuda_implementation/cuda_random.cuh>
#include <chrono>
#include <mudock/cuda_implementation/cuda_check_error_macro.cuh>

namespace mudock {
  __global__ void init_curand(curandState *state, const long seed, const int num_elements) {
    const int id     = threadIdx.x + blockIdx.x * blockDim.x;
    const int stride = gridDim.x * blockDim.x;
    for (int index = id; index < num_elements; index += stride) {
      curand_init(seed + index, index, 0, &state[index]);
    }
  }

  void cuda_random_object::alloc(const std::size_t num_elements) {
    const bool init = num_elements > cuda_object<curandState>::num_elements();
    cuda_object<curandState>::alloc(num_elements);
    if (init) {
      init_curand<<<4, 128>>>(cuda_object<curandState>::dev_pointer(),
                              std::chrono::high_resolution_clock::now().time_since_epoch().count(),
                              cuda_object<curandState>::num_elements());
      MUDOCK_CHECK_KERNELCALL();
      MUDOCK_CHECK(cudaDeviceSynchronize());
    }
  };
} // namespace mudock
