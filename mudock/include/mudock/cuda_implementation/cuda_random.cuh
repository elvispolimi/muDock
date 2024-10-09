#pragma once

#include <mudock/cuda_implementation/cuda_wrapper.cuh>
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

    // Generate the next random number
    __host__ __device__ fp_type next() {
      /* Algorithm "xorwow" from p. 5 of Marsaglia, "Xorshift RNGs" */
      unsigned int t = state[4];

      const unsigned int s = state[0]; /* Perform a contrived 32-bit rotate. */
      state[4]             = state[3];
      state[3]             = state[2];
      state[2]             = state[1];
      state[1]             = s;

      t ^= t >> 2;
      t ^= t << 1;
      t ^= s ^ (s << 4);
      state[0] = t;
      index += 362437;
      return static_cast<fp_type>(t + index) / static_cast<fp_type>(std::numeric_limits<unsigned int>::max());
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
