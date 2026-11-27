#pragma once

#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/compute/geom_transform_args.hpp>
#include <mudock/type_alias.hpp>
#include <cstddef>
#include <type_traits>

namespace mudock {
  struct genetic_args {
    int batch_ligands;
    int population_size;
    int generations;
    std::size_t seed;
    int tournament_length;
    fp_type mutation_probability;
    fp_type max_translation_half_width;
    chromosome* population;
    chromosome* next_population;
    const int* num_rotamers;
    fp_type* scores;
    fp_type* best_scores;
    chromosome* best_chromosomes;
    geom_transform_args geometry;
  };
  static_assert(std::is_trivially_copyable_v<genetic_args>);
} // namespace mudock
