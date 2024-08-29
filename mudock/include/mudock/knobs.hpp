#pragma once

#include <cstdint>

namespace mudock {

  struct knobs {
    std::size_t population_number = 100;
    std::size_t num_generations   = 1000;
    std::size_t tournament_length = 10;
    std::size_t mutation_prob     = 0.01;
  };

} // namespace mudock
