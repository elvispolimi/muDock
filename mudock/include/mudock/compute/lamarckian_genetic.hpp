#pragma once

#include "mudock/knobs.hpp"
#include "genetic.hpp"

#include <concepts>
#include <memory>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_protein.hpp>
#ifndef __CUDACC__
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/docking.hpp>
  #include <mudock/compute/geometric_transform.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/compute/queue.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/molecule.hpp>  

namespace mudock {

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct lamarckian_genetic_kernel : public genetic_kernel<queue_type> {
    static constexpr char finalize_region_name[]   = "lamarckian_genetic_finalize";
    static constexpr char iterate_region_name[]    = "lamarckian_genetic_iterate";
    static constexpr char initialize_region_name[] = "lamarckian_genetic_initialize";

    lamarckian_genetic_kernel(const int batch_ligands_,
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
                              chromosome* __restrict__ best_chromosomes_b_,
                              std::shared_ptr<queue_type> q_)
        : genetic_kernel<queue_type>(batch_ligands_,
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
    lamarckian_genetic_kernel() {};

    void operator()() { genetic_kernel<queue_type>::operator()(); }
    void initialize() { genetic_kernel<queue_type>::initialize(); }
    void finalize() { genetic_kernel<queue_type>::finalize(); }

  private:
  };

#ifndef __CUDACC__
  template<typename queue_t, template<typename> typename scoring_t, template<typename> typename local_search_t>
    requires std::derived_from<queue_t, queue> 
             && std::derived_from<scoring_t<queue_t>, scoring<queue_t>>
             && std::derived_from<local_search_t<queue_t>, local_search<queue_t>>
  struct lamarckian_genetic : public genetic<queue_t, scoring_t> {
    // Inherit constructor and most logic from genetic
    lamarckian_genetic(std::shared_ptr<scratchpad<queue_t>> _scratch,
                      dynamic_molecule& _protein,
                      scoring_t<queue_t> _scoring,
                      std::unique_ptr<local_search> _local_search,
                      int _local_search_iters)
        : genetic<queue_t, scoring_t>(_scratch, _protein, _scoring),
          local_search_stage(std::move(_local_search)),
          local_search_iterations(_local_search_iters) {}

    void prepare(batch<static_molecule>& batch) {
      const knobs& configuration  = (*this->scratch).configuration;
      batch_ligands               = batch.num_ligands;
      num_generations             = configuration.num_generations;
      const int population_number = configuration.population_number;
      auto q                      = (*this->scratch).get_queue();

      auto& num_rotamers_b = (*this->scratch).template get<buffer_data_type::NUM_ROTAMERS>();
      auto& chromosomes_b  = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      auto& scores_b       = (*this->scratch).template get<buffer_data_type::SCORES>();

      num_rotamers_b.alloc(batch_ligands);
      chromosomes_b.alloc(population_number * batch_ligands);
      next_population.alloc(population_number * batch_ligands);
      scores_b.alloc(population_number * batch_ligands);
      best_scores.alloc(batch_ligands);
      best_chromosomes.alloc(batch_ligands);

      load_num_rotamers<queue_t>(batch, this->scratch);

      const auto seed =
          configuration.seed.has_value()
              ? configuration.seed.value()
              : static_cast<size_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());

      int* __restrict__ num_rotamers_p            = num_rotamers_b.dev_pointer();
      fp_type* __restrict__ scores_p              = scores_b.dev_pointer();
      fp_type* __restrict__ best_scores_p         = best_scores.dev_pointer();
      chromosome* __restrict__ best_chromosomes_p = best_chromosomes.dev_pointer();

      this->lamarckian_kernel = std::make_unique<lamarckian_genetic_kernel<queue_t>>(batch_ligands,
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
                                                                                     best_chromosomes_p,
                                                                                     q);

      geom_trans.prepare(batch);
      local_search_stage.prepare(batch);
      score_stage.prepare(batch);
    };

    void operator()() override {
      // TODO try to recycle genetic code and adding just local search 
      auto& chromosomes_b = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();

      assert(lamarckian_kernel && "lamarckian_kernel method not yet prepared");
      lamarckian_kernel->initialize();

      for (int generation = 0; generation < num_generations; ++generation) {
        geom_trans();
        local_search_stage();
        score_stage();
        (*lamarckian_kernel)();
        chromosomes_b.copy_device2device(next_population);
      }
      lamarckian_kernel->finalize();
    }

  private:
    int local_search_iterations;
    local_search_t<queue_t> local_search_stage;
    std::unique_ptr<lamarckian_genetic_kernel<queue_t>> lamarckian_kernel;
    
    // TODO L probably this must be changed, not sure
    void teardown_impl(batch<static_molecule>& batch) {
      assert(batch.num_ligands == batch_ligands && "Lamarckian-Genetic algorithm received different batch for teardown");

      // TODO check if the copy can be changed with a swap
      // auto& population_b = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      // population_b.copy_device2device(best_chromosomes);
      // geom_trans();
      // score_stage();
      // geom_trans.teardown(batch);
      // score_stage.teardown(batch);

      best_scores.copy_device2host();
      (*this->scratch).get_queue()->synchronize();

      for (int index{0}; index < batch_ligands; ++index) {
        auto& ligand = *batch.molecules[index];
        ligand.properties.assign(property_type::SCORE, std::to_string(best_scores()[index]));
      }
    }
  }; // namespace mudock
#endif
} // namespace mudock
