#include <cstdint>
#include <mudock/hip_implementation/hip_random.hpp>
#include <mudock/hip_implementation/hip_utils.hpp>

namespace mudock {
  namespace {
    std::size_t default_hip_seed() {
      // Keep HIP random seeding independent of <chrono> to avoid dragging the
      // GCC 13 chrono/format header chain into HIP compilation.
      static std::uint64_t counter = 0x9e3779b97f4a7c15ULL;
      counter += 0x9e3779b97f4a7c15ULL;
      return static_cast<std::size_t>(counter);
    }
  } // namespace

  __global__ void init_hiprand(hiprandState *state, const long seed, const int num_elements) {
    const int id     = threadIdx.x + blockIdx.x * blockDim.x;
    const int stride = gridDim.x * blockDim.x;
    for (int index = id; index < num_elements; index += stride) {
      hiprand_init(seed + index, index, 0, &state[index]);
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
    alloc(num_elements, default_hip_seed());
  };
} // namespace mudock
