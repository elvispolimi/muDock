#include <algorithm>
#include <chrono>
#include <cstddef>
#include <mudock/alpaka_implementation/alpaka_random.hpp>
#include <mudock/alpaka_implementation/invoke_kernel_alpaka.hpp>

namespace mudock {

  struct init_alpaka_rand {
    template<typename TAcc>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  alpaka_rand_state* state,
                                  const std::size_t seed,
                                  const std::size_t num_elements) const {
      const int id     = static_cast<int>(alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u]);
      const int stride = static_cast<int>(alpaka::getWorkDiv<alpaka::Grid, alpaka::Threads>(acc)[0u]);

      for (std::size_t i = id; i < num_elements; i += stride) { state[i] = alpaka_rand_state(seed, i, 0); }
    }
  };

  void alpaka_random_object::alloc(const std::size_t num_elements, const std::size_t seed) {
    const auto current_size = state.num_elements();
    const bool needs_init   = num_elements > current_size;

    state.alloc(num_elements);

    if (needs_init) {
      q->invoke_kernel<init_alpaka_rand>(128, 32, state.dev_pointer(), seed, num_elements);
      q->synchronize();
    }
  }

  void alpaka_random_object::alloc(const std::size_t num_elements) {
    alloc(num_elements, std::chrono::high_resolution_clock::now().time_since_epoch().count());
  }

} // namespace mudock
