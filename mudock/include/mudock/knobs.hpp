#pragma once

#include <cstddef>
#include <mudock/type_alias.hpp>
#include <optional>

namespace mudock {

  struct knobs {
    std::size_t population_number   = 100;
    std::size_t num_generations     = 1000;
    std::size_t tournament_length   = 10;
    fp_type mutation_prob           = static_cast<fp_type>(0.01);
    std::optional<std::size_t> seed = std::optional<std::size_t>{};
    std::size_t max_tbb_tokens      = 4;
    std::size_t max_bytes_per_token = 1048576; // 1 MB
    std::size_t max_tbb_queue_size  = 100000;
  };

} // namespace mudock
