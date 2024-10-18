#pragma once

#include <polygeist/cuda_wrapper.cuh>
#include <mudock/type_alias.hpp>
#include <random>

namespace mudock {
  // XORWOW state
  // Should be the default one from curand state in CUDA
  // Taken from Wikipedia https://en.wikipedia.org/wiki/Xorshift#xorwow
  struct XORWOWState {
    std::array<unsigned int, 5> state; // State array
    int index;                         // Current index in the state array

    XORWOWState() = default;

    // Constructor to initialize the state with a seed
    void set_seed(const int seed) {
      index    = 5;
      state[0] = seed;
      for (unsigned int i = 1; i < 5; ++i) {
        state[i] = 1812433253U * (state[i - 1] ^ (state[i - 1] >> 30)) + i;
      }
    }
  };

  struct cuda_random_object: private cuda_wrapper<std::vector, XORWOWState> {
    cuda_random_object()                                      = default;
    cuda_random_object(const cuda_random_object &)            = delete;
    cuda_random_object(cuda_random_object &&)                 = default;
    cuda_random_object &operator=(const cuda_random_object &) = delete;
    cuda_random_object &operator=(cuda_random_object &&)      = default;

    void alloc(const std::size_t num_elements) {
      const bool init = num_elements > cuda_wrapper<std::vector, XORWOWState>::num_elements();
      cuda_wrapper<std::vector, XORWOWState>::alloc(num_elements);
      auto generator = std::mt19937{
          static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count())};
      auto dist = std::uniform_int_distribution();
      if (init) {
        for (auto &s: host) { s.set_seed(dist(generator)); }
        cuda_wrapper<std::vector, XORWOWState>::copy_host2device();
      }
    };

    [[nodiscard]] inline auto dev_pointer() const {
      return cuda_wrapper<std::vector, XORWOWState>::dev_pointer();
    }
  };

} // namespace mudock
