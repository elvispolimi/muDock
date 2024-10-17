#include <hiprand/hiprand_kernel.h>
#include <mudock/hip_implementation/hip_random.hpp>

namespace mudock {
  __global__ void init_hiprand(hiprandState *state, const long seed, const int num_elements) {
    const int id     = threadIdx.x + blockIdx.x * blockDim.x;
    const int stride = gridDim.x * blockDim.x;
    for (int index = id; index < num_elements; index += stride) {
      // const int seed = hash(seed+index); // Create a unique seed for each thread
      hiprand_init(seed + index, index, 0, &state[index]);
    }
  }

  void hip_random_object::alloc(const std::size_t num_elements) {
    const bool init = num_elements > hip_object<hiprandState>::num_elements();
    hip_object<hiprandState>::alloc(num_elements);
    if (init) {
      init_hiprand<<<4, 128, 0, hip_object<hiprandState>::get_stream()>>>(
          hip_object<hiprandState>::dev_pointer(),
          std::chrono::high_resolution_clock::now().time_since_epoch().count(),
          hip_object<hiprandState>::num_elements());
      MUDOCK_CHECK_KERNELCALL();
      MUDOCK_CHECK(hipStreamSynchronize(hip_object<hiprandState>::get_stream()));
    }
  };
} // namespace mudock
