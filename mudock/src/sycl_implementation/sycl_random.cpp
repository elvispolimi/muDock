#include <mudock/sycl_implementation/sycl_random.hpp>
#include <random>

namespace mudock {

  void sycl_random_object::alloc(const std::size_t num_elements, const std::size_t seed) {
    const bool init = num_elements > sycl_wrapper<std::vector, XORWOWState>::num_elements();
    sycl_wrapper<std::vector, XORWOWState>::alloc(num_elements);
    auto generator = std::mt19937{static_cast<unsigned>(seed)};
    auto dist      = std::uniform_int_distribution();
    if (init) {
      for (auto& s: host) { s.set_seed(dist(generator)); }
      sycl_wrapper<std::vector, XORWOWState>::copy_host2device();
    }
  };
  void sycl_random_object::alloc(const std::size_t num_elements) {
    alloc(num_elements, std::chrono::high_resolution_clock::now().time_since_epoch().count());
  };
} // namespace mudock
