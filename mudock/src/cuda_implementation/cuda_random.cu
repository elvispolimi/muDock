#include <mudock/cuda_implementation/cuda_random.cuh>

namespace mudock {
  __global__ void init_curand(curandState *state, const int num_elements) {
    const int id     = threadIdx.x + blockIdx.x * blockDim.x;
    const int stride = gridDim.x * blockDim.x;
    for (int index = id; index < num_elements; index += stride) curand_init(clock64(), id, 0, &state[id]);
  }

  void cuda_random_object::alloc(const std::size_t num_elements) {
    const bool init = num_elements > cuda_object<curandState>::num_elements();
    cuda_object<curandState>::alloc(num_elements);
    if (init) {
      init_curand<<<32, 128>>>(cuda_object<curandState>::dev_pointer(),
                               cuda_object<curandState>::num_elements());
      MUDOCK_CHECK_KERNELCALL();
      MUDOCK_CHECK(cudaDeviceSynchronize());
    }
  };
} // namespace mudock