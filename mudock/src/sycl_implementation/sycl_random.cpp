#include <chrono>
#include <mudock/sycl_implementation/sycl_random.hpp>
#include <random>

namespace mudock {

  void sycl_random_object::alloc(const std::size_t num_elements, const std::size_t seed) {
    const bool init = num_elements > state.num_elements();
    state.alloc(num_elements);
    if (init) {
      auto generator = std::mt19937{static_cast<unsigned>(seed)};
      auto dist      = std::uniform_int_distribution();
      for (size_t i = 0; i < num_elements; i++) state()[i].set_seed(dist(generator));
      state.copy_host2device();
      q->synchronize();
    }
  };
  void sycl_random_object::alloc(const std::size_t num_elements) {
    alloc(num_elements, std::chrono::high_resolution_clock::now().time_since_epoch().count());
  };
} // namespace mudock
