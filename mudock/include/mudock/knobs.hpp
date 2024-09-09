#pragma once

#include "mudock/type_alias.hpp"

#include <cstdint>

namespace mudock {

  struct knobs {
    std::size_t population_number = 100;
    std::size_t num_generations   = 1000;
    std::size_t tournament_length = 10;
    fp_type mutation_prob         = fp_type{0.01};
  };

} // namespace mudock
