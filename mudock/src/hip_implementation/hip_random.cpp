#include <chrono>
#include <mudock/hip_implementation/hip_random.hpp>
#include <mudock/hip_implementation/hip_utils.hpp>

namespace mudock {
  __global__ void init_hiprand(hiprandState *state, const long seed, const int num_elements) {
    const int id     = threadIdx.x + blockIdx.x * blockDim.x;
    const int stride = gridDim.x * blockDim.x;
    for (int index = id; index < num_elements; index += stride) {
      curand_init(seed + index, index, 0, &state[index]);
    }
  }

  void hip_random_object::alloc(const std::size_t num_elements, const std::size_t seed) {
    const auto num_el = state.num_elements();
    //TODO check this condition
    const bool init = num_elements > num_el;
    state.alloc(num_elements);
    if (init) {
      void *args[] = {(void *) state.dev_pointer_ref(), (void *) &seed, (void *) &num_elements};
      //TODO check grid/dimensions
      q->launch_kernel((void *) init_hiprand, args, 128, BLOCK_SIZE);
      q->synchronize();
    }
  };

  void hip_random_object::alloc(const std::size_t num_elements) {
    alloc(num_elements, std::chrono::high_resolution_clock::now().time_since_epoch().count());
  };
} // namespace mudock
