#include <mudock/sycl_implementation/sycl_random.hpp>

namespace mudock {
  unsigned int hash(unsigned int x) {
    return (x ^ (x >> 16)) * 0x45d9f301; // A simple hash function
  }

  float uniform_real(XORWOWState& state) {
    return static_cast<float>(state.next()) / static_cast<float>(0xFFFFFFFFU);
  }

  void sycl_random_object::alloc(const std::size_t num_elements) {
    const bool init = num_elements > sycl_wrapper<std::vector, XORWOWState>::num_elements();
    sycl_wrapper<std::vector, XORWOWState>::alloc(num_elements);
    // TODO maybe initialize brand new
    if (init) {
      int index{0};
      for (auto& s: host) {
        s.set_seed(hash(index));
        ++index;
      }
      sycl_wrapper<std::vector, XORWOWState>::copy_host2device();
    }
  };
} // namespace mudock
