#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/genetic.hpp>
#include <mudock/cpp_implementation/genetic_cpp.hpp>
#include <mudock/gh_implementation/queue_gh.hpp>

namespace mudock {

  template<>
  struct genetic_kernel<queue_gh>: genetic_kernel<queue_cpp> {
    genetic_kernel(const int batch_ligands_,
                   const int population_number_,
                   const int num_generations_,
                   const int tournament_length_,
                   const fp_type mutation_prob_,
                   const size_t seed_,
                   chromosome* population_,
                   chromosome* next_population_,
                   int* __restrict__ num_rotamers_b_,
                   fp_type* __restrict__ scores_b_,
                   fp_type* __restrict__ best_scores_b_,
                   chromosome* __restrict__ best_chromosomes_b_,
                   std::shared_ptr<queue_gh> q_)
        : genetic_kernel<queue_cpp>(batch_ligands_,
                                    population_number_,
                                    num_generations_,
                                    tournament_length_,
                                    mutation_prob_,
                                    seed_,
                                    population_,
                                    next_population_,
                                    num_rotamers_b_,
                                    scores_b_,
                                    best_scores_b_,
                                    best_chromosomes_b_,
                                    q_) {};
    genetic_kernel() {};
  };
} // namespace mudock
