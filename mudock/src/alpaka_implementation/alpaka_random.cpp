/**
 * @file alpaka_random.cpp
 * @brief Implementation of device-side PRNG initialization and allocation.
 * @details Implements the `init_alpaka_rand` initialization kernel and methods of
 *          `alpaka_random_object` to seed Philox/Mersenne-Twister engines on accelerator devices.
 */

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <mudock/alpaka_implementation/alpaka_random.hpp>
#include <mudock/alpaka_implementation/invoke_kernel_alpaka.hpp>

namespace mudock {

  /**
   * @struct init_alpaka_rand
   * @brief Kernel initializing device PRNG states across a 1D grid.
   * @details Computes global thread ID and grid stride to initialize `alpaka_rand_state`
   *          instances with unique sub-sequences derived from a common base seed and sequence index.
   */
  struct init_alpaka_rand {
    /**
     * @brief Kernel body initializing random states.
     * @tparam TAcc Alpaka accelerator type.
     * @param[in] acc Reference to the execution context.
     * @param[out] state Device array of random generator states.
     * @param seed Base numeric seed.
     * @param num_elements Total number of generator states to initialize.
     */
    template<typename TAcc>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  alpaka_rand_state* state,
                                  const std::size_t seed,
                                  const std::size_t num_elements) const {
      const int id     = static_cast<int>(alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u]);
      const int stride = static_cast<int>(alpaka::getWorkDiv<alpaka::Grid, alpaka::Threads>(acc)[0u]);

      for (std::size_t i = id; i < num_elements; i += stride) { 
        state[i] = alpaka_rand_state(static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(i), 0); 
      }
    }
  };

  /**
   * @brief Allocates and initializes device random state buffers with an explicit seed.
   * @param num_elements Required number of generator states.
   * @param seed Seed value used to initialize states if buffer is resized.
   */
  void alpaka_random_object::alloc(const std::size_t num_elements, const std::size_t seed) {
    const auto current_size = state.num_elements();
    const bool needs_init   = num_elements > current_size;

    state.alloc(num_elements);

    if (needs_init) {
      q->invoke_kernel<init_alpaka_rand>(128, state.dev_pointer(), seed, num_elements);
      q->synchronize();
    }
  }

  /**
   * @brief Allocates and initializes device random state buffers using high-resolution clock seed.
   * @param num_elements Required number of generator states.
   */
  void alpaka_random_object::alloc(const std::size_t num_elements) {
    alloc(num_elements, std::chrono::high_resolution_clock::now().time_since_epoch().count());
  }

} // namespace mudock
