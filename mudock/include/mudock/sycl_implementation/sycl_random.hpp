#pragma once

#include <mudock/type_alias.hpp>
#include <array>
#include <limits>

namespace mudock {
  // Lightweight per-work-item XORWOW state for SYCL.
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
    inline fp_type next() {
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

  // Use vendor-style per-work-item seed/sequence initialization and sampling.
  // The custom XORWOW state folds sequence into the seed; it has no jump-ahead
  // offset like CURAND/HIPRAND.
  inline void syclrand_init(const unsigned long long seed,
                            const unsigned long long sequence,
                            XORWOWState* state) {
    state->set_seed(static_cast<int>(seed + sequence));
  }

  inline fp_type syclrand_uniform(XORWOWState* state) { return state->next(); }

} // namespace mudock
