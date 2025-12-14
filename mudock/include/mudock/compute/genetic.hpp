#pragma once

#include <concepts>
#include <memory>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/buffer.hpp>
#include <mudock/compute/docking.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/compute/scoring.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct rand_state_type;

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct genetic_kernel {
    using rand_type = typename rand_state_type<queue_type>::type;

    genetic_kernel(const int batch_ligands_,
                   const int population_number_,
                   const int num_generations_,
                   const int tournament_length_,
                   const int mutation_prob_,
                   const size_t seed_,
                   chromosome* population_,
                   chromosome* next_population_,
                   int* __restrict__ num_rotamers_b_,
                   fp_type* __restrict__ scores_b_,
                   fp_type* __restrict__ best_scores_b_,
                   chromosome* __restrict__ best_chromosomes_b_)
        : batch_ligands(batch_ligands_),
          population_number(population_number_),
          num_generations(num_generations_),
          tournament_length(tournament_length_),
          mutation_prob(mutation_prob_),
          population(population_),
          next_population(next_population_),
          num_rotamers_b(num_rotamers_b_),
          scores_b(scores_b_),
          best_scores_b(best_scores_b_),
          best_chromosomes_b(best_chromosomes_b_),
          rand(seed_) {};
    genetic_kernel() {};

    void operator()();
    void initialize();
    void finalize();

  private:
    int batch_ligands;
    int population_number;
    int num_generations;
    int tournament_length;
    int mutation_prob;
    chromosome* __restrict__ population;
    chromosome* __restrict__ next_population;
    int* __restrict__ num_rotamers_b;
    fp_type* __restrict__ scores_b;
    fp_type* __restrict__ best_scores_b;
    chromosome* __restrict__ best_chromosomes_b;
    rand_type rand;
  };

  template<typename queue_t, template<typename> typename scoring_t>
    requires std::derived_from<queue_t, queue> && std::derived_from<scoring_t<queue_t>, scoring<queue_t>>
  struct genetic: public docking<queue_t> {
    genetic(std::shared_ptr<scratchpad<queue_t>> _scratch,
            dynamic_molecule& protein,
            scoring_t<queue_t> _scoring)
        : docking<queue_t>(_scratch),
          score_stage(std::move(_scoring)),
          geom_trans(_scratch, protein),
          next_population(_scratch->get_queue()),
          best_chromosomes(_scratch->get_queue()),
          best_scores(_scratch->get_queue()) {};
    void prepare(batch<static_molecule>& batch) {
      const knobs& configuration  = (*this->scratch).configuration;
      batch_ligands               = batch.num_ligands;
      num_generations             = configuration.num_generations;
      const int population_number = configuration.population_number;

      auto& num_rotamers_b = (*this->scratch).template get<buffer_data_type::NUM_ROTAMERS>();
      auto& chromosomes_b  = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      auto& scores_b       = (*this->scratch).template get<buffer_data_type::SCORES>();

      num_rotamers_b.alloc(batch_ligands);
      chromosomes_b.alloc(population_number * batch_ligands);
      next_population.alloc(population_number * batch_ligands);
      scores_b.alloc(population_number * batch_ligands);
      best_scores.alloc(batch_ligands);
      best_chromosomes.alloc(batch_ligands);

      for (int index{0}; index < batch_ligands; ++index) {
        auto& ligand            = *batch.molecules[0];
        num_rotamers_b()[index] = ligand.num_rotamers();
      }

      num_rotamers_b.copy_host2device();

      const auto seed =
          configuration.seed.has_value()
              ? configuration.seed.value()
              : static_cast<size_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());

      int* __restrict__ num_rotamers_p            = num_rotamers_b.dev_pointer();
      fp_type* __restrict__ scores_p              = scores_b.host_pointer();
      fp_type* __restrict__ best_scores_p         = best_scores.dev_pointer();
      chromosome* __restrict__ best_chromosomes_p = best_chromosomes.dev_pointer();

      kernel = std::make_unique<genetic_kernel<queue_t>>(batch_ligands,
                                                         population_number,
                                                         configuration.num_generations,
                                                         configuration.tournament_length,
                                                         configuration.mutation_prob,
                                                         seed,
                                                         chromosomes_b.dev_pointer(),
                                                         next_population.dev_pointer(),
                                                         num_rotamers_p,
                                                         scores_p,
                                                         best_scores_p,
                                                         best_chromosomes_p);

      geom_trans.prepare(batch);
      score_stage.prepare(batch);
    };
    void operator()() {
      auto& chromosomes_b = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      kernel->initialize();

      for (int generation = 0; generation < num_generations; ++generation) {
        geom_trans();
        score_stage();
        (*kernel)();
        chromosomes_b.copy_device2device(next_population);
      }
      kernel->finalize();
    };

  private:
    // std::shared_ptr<scoring<queue_t>> score_stage;
    scoring_t<queue_t> score_stage;
    geometric<queue_t> geom_trans;
    std::unique_ptr<genetic_kernel<queue_t>> kernel;

    int batch_ligands;
    int num_generations;
    buffer_vector<chromosome, queue_t> next_population;
    buffer_vector<chromosome, queue_t> best_chromosomes;
    buffer_vector<fp_type, queue_t> best_scores;

    void teardown_impl(batch<static_molecule>& batch) {
      assert(batch.num_ligands == batch_ligands && "Genetic algorithm received different batch for teardown");

      auto& population_b       = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      auto& best_chromosomes_b = best_chromosomes;

      // TODO check if the copy can be changed with a swap
      // population_b.swap(best_chromosomes_b);
      population_b.copy_device2device(best_chromosomes_b);

      geom_trans();
      // geom_trans.teardown(batch);
      score_stage();
      score_stage.teardown(batch);

      for (int index{0}; index < batch_ligands; ++index) {
        auto& ligand = *batch.molecules[index];
        ligand.properties.assign(property_type::SCORE, std::to_string(best_scores()[index]));
      }
    }
  }; // namespace mudock
} // namespace mudock
