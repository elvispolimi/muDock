#pragma once

#include "mudock/type_alias.hpp"

#include <cstdint>

namespace mudock {

  struct knobs {
    std::size_t population_number = 10;
    std::size_t num_generations   = 100;
    std::size_t tournament_length = 10;
    fp_type mutation_prob         = fp_type{0.01};
  };

} // namespace mudock
