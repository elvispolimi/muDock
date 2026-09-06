#pragma once

#include <cstddef>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <optional>

namespace mudock {

  struct knobs {
    std::size_t population_number             = 100;
    std::size_t num_generations               = 1000;
    std::size_t elite_size                    = 0;
    std::size_t tournament_length             = 10;
    fp_type mutation_prob                     = static_cast<fp_type>(0.01);
    bool autostop                             = false;
    std::size_t tolerance_window              = 50;
    fp_type score_variance_thld               = static_cast<fp_type>(0.0015);
    fp_type best_score_diff_thld              = static_cast<fp_type>(0.001);
    // fp_type crystal_score                     = small_bound<fp_type>();
    // fp_type crystal_tolerance                 = static_cast<fp_type>(0.5);
    fp_type lsrate                            = static_cast<fp_type>(100);
    bool ls_on_best                           = false;
    std::size_t lsit                          = 300;
    std::size_t ls_every                      = 1;
    std::size_t ls_last_gen                   = 0;
    std::optional<std::size_t> seed           = std::optional<std::size_t>{};
    std::size_t max_tbb_tokens                = 4;
    std::size_t max_bytes_per_token           = 1048576; // 1 MB
    std::size_t max_tbb_queue_size            = 100000;
  };

} // namespace mudock
