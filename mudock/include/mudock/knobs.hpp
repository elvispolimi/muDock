#pragma once

#include <cstddef>
#include <mudock/type_alias.hpp>
#include <optional>

namespace mudock {

  struct knobs {
    std::size_t population_number                      = 100;
    std::size_t num_generations                        = 1000;
    static constexpr std::size_t num_score_generations = 1;
    std::size_t tournament_length                      = 10;
    fp_type mutation_prob                              = fp_type{0.01};
    std::optional<std::size_t> seed                    = std::optional<std::size_t>{};
  };

} // namespace mudock
