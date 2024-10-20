#pragma once

#include <array>
#include <mudock/omp_implementation/omp_wrapper.hpp>
#include <mudock/type_alias.hpp>
#include <vector>

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

  struct omp_random_object: private omp_wrapper<std::vector, XORWOWState> {
    omp_random_object(): omp_wrapper<std::vector, XORWOWState>(){};
    omp_random_object(const omp_random_object &)            = delete;
    omp_random_object(omp_random_object &&)                 = default;
    omp_random_object &operator=(const omp_random_object &) = delete;
    omp_random_object &operator=(omp_random_object &&)      = delete;

    void alloc(const std::size_t num_elements);

    inline void copy_host2device() { omp_wrapper<std::vector, XORWOWState>::copy_host2device(); }

    [[nodiscard]] inline auto dev_pointer() const {
      return omp_wrapper<std::vector, XORWOWState>::dev_pointer();
    }
  };

} // namespace mudock
