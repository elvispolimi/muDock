#pragma once

#include "mudock/knobs.hpp"

#include <algorithm>
#include <concepts>
#include <limits>
#include <memory>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <numeric>
#include <stdexcept>
#if !defined(__CUDACC__) && !defined(__HIPCC__)
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/docking.hpp>
  #include <mudock/compute/geometric_transform.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct genetic_kernel {
    static constexpr char finalize_region_name[]   = "genetic_finalize";
    static constexpr char iterate_region_name[]    = "genetic_iterate";
    static constexpr char initialize_region_name[] = "genetic_initialize";

    genetic_kernel(const int batch_ligands_,
                   const int population_number_,
                   const int num_generations_,
                   const int elite_size_,
                   const fp_type score_variance_thld_,
                   const fp_type best_score_diff_thld_,
                  //  const fp_type crystal_score_,
                  //  const fp_type crystal_tolerance_,
                   const bool autostop_,
                   int* __restrict__ converged_ligands_b_,
                   fp_type* __restrict__ best_so_far_b_,
                   int* __restrict__ for_how_long_best_b_,
                   const int tolerance_window_,
                   const int tournament_length_,
                   const fp_type mutation_prob_,
                   const size_t seed_,
                   chromosome* population_,
                   chromosome* next_population_,
                   int* __restrict__ num_rotamers_b_,
                   fp_type* __restrict__ scores_b_,
                   fp_type* __restrict__ best_scores_b_,
                   chromosome* __restrict__ best_chromosomes_b_,
                   std::shared_ptr<queue_type> q_)
        : batch_ligands(batch_ligands_),
          population_number(population_number_),
          num_generations(num_generations_),
          elite_size(elite_size_),
          score_variance_thld(score_variance_thld_),
          best_score_diff_thld(best_score_diff_thld_),
          // crystal_score(crystal_score_),
          // crystal_tolerance(crystal_tolerance_),
          autostop(autostop_),
          converged_ligands_b(converged_ligands_b_),
          best_so_far_b(best_so_far_b_),
          for_how_long_best_b(for_how_long_best_b_),
          tolerance_window(tolerance_window_),
          tournament_length(tournament_length_),
          mutation_prob(mutation_prob_),
          population(population_),
          next_population(next_population_),
          num_rotamers_b(num_rotamers_b_),
          scores_b(scores_b_),
          best_scores_b(best_scores_b_),
          best_chromosomes_b(best_chromosomes_b_),
          seed(seed_),
          q(q_) {};
    genetic_kernel() {};

    void operator()();
    void initialize();
    void finalize();
    inline void set_population_buffers(chromosome* population_, chromosome* next_population_) {
      population      = population_;
      next_population = next_population_;
    }

  private:
    int batch_ligands;
    int population_number;
    int num_generations;
    int elite_size;
    fp_type score_variance_thld;  
    fp_type best_score_diff_thld;  
    // fp_type crystal_score;
    // fp_type crystal_tolerance;
    bool autostop;
    int* __restrict__ converged_ligands_b;
    fp_type* __restrict__ best_so_far_b;
    int* __restrict__ for_how_long_best_b;
    int tolerance_window;
    int tournament_length;
    int current_generation = 1;
    fp_type mutation_prob;
    chromosome* __restrict__ population;
    chromosome* __restrict__ next_population;
    int* __restrict__ num_rotamers_b;
    fp_type* __restrict__ scores_b;
    fp_type* __restrict__ best_scores_b;
    chromosome* __restrict__ best_chromosomes_b;
    size_t seed;
    std::shared_ptr<queue_type> q;
  };

#if !defined(__CUDACC__) && !defined(__HIPCC__)
  template<typename queue_t, template<typename> typename scoring_t>
    requires std::derived_from<queue_t, queue> && std::derived_from<scoring_t<queue_t>, scoring<queue_t>>
  struct genetic: public docking<queue_t> {
    static constexpr const char stage_name[] = "GENETIC";

    genetic(std::shared_ptr<scratchpad<queue_t>> _scratch,
            dynamic_molecule& protein,
            std::shared_ptr<scoring_t<queue_t>> _scoring)
        : docking<queue_t>(_scratch),
          score_stage(std::move(_scoring)),
          geom_trans(_scratch, protein),
          next_population(_scratch->get_queue()),
          best_chromosomes(_scratch->get_queue()),
          best_scores(_scratch->get_queue()) {};
    void prepare(batch<static_molecule>& batch) {
      const knobs& configuration   = (*this->scratch).configuration;
      batch_ligands                = batch.num_ligands;
      num_generations              = static_cast<int>(configuration.num_generations);
      population_number            = static_cast<int>(configuration.population_number);
      auto q                       = (*this->scratch).get_queue();

      auto& num_rotamers_b            = (*this->scratch).template get<buffer_data_type::NUM_ROTAMERS>();
      auto& chromosomes_b             = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      auto& scores_b                  = (*this->scratch).template get<buffer_data_type::SCORES>();
      auto& converged_ligands_b       = (*this->scratch).template get<buffer_data_type::CONVERGED_LIGANDS>();
      auto& best_so_far_b             = (*this->scratch).template get<buffer_data_type::BEST_SO_FAR>();
      auto& for_how_long_best_b       = (*this->scratch).template get<buffer_data_type::FOR_HOW_LONG_BEST>();

      num_rotamers_b.alloc(batch_ligands);
      chromosomes_b.alloc(population_number * batch_ligands);
      next_population.alloc(population_number * batch_ligands);
      scores_b.alloc(population_number * batch_ligands);
      best_scores.alloc(batch_ligands);
      best_chromosomes.alloc(batch_ligands);
      converged_ligands_b.alloc(batch_ligands);
      best_so_far_b.alloc(batch_ligands);
      for_how_long_best_b.alloc(batch_ligands);

      load_num_rotamers<queue_t>(batch, this->scratch);
      initialize_converged_ligands<queue_t>(batch, this->scratch);
      initialize_best_so_far<queue_t>(batch, this->scratch);
      initialize_for_how_long_best<queue_t>(batch, this->scratch);

      const auto seed =
          configuration.seed.has_value()
              ? configuration.seed.value()
              : static_cast<size_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());

      int* __restrict__ num_rotamers_p                = num_rotamers_b.dev_pointer();
      fp_type* __restrict__ scores_p                  = scores_b.dev_pointer();
      fp_type* __restrict__ best_scores_p             = best_scores.dev_pointer();
      chromosome* __restrict__ best_chromosomes_p     = best_chromosomes.dev_pointer();
      int* __restrict__ converged_ligands_p           = converged_ligands_b.dev_pointer();
      fp_type* __restrict__ best_so_far_p             = best_so_far_b.dev_pointer();
      int* __restrict__ for_how_long_best_p           = for_how_long_best_b.dev_pointer();
      

      // Lorenzo: Ligand properties for experiments
      for (int index{0}; index < batch_ligands; ++index) {
        auto& ligand = *batch.molecules[index];
        const int num_rotamers = num_rotamers_p[index];
        // TODO L check if this estimate is correct. WARNING: this depends on the local search implementation. 
        // This is for adadelta for example (not counting the effect of ls_last_gen and ls_every)
        // +1 comes from the scores, the remaining from the gradients
        // TODO L this is not correct: if autostop is on, it should count the actual number of generations at convergence.
        // It would be better to have a counter at each evaluation to be sure (pay attention to race conditions)
        ligand.properties.assign(property_type::POP, std::to_string(population_number));
        ligand.properties.assign(property_type::SEED, std::to_string(seed));
        ligand.properties.assign(property_type::NUM_ROT, std::to_string(num_rotamers));
      }

      kernel = std::make_unique<genetic_kernel<queue_t>>(batch_ligands,
                                                         population_number,
                                                         configuration.num_generations,
                                                         configuration.elite_size,
                                                         configuration.score_variance_thld,
                                                         configuration.best_score_diff_thld,
                                                        //  configuration.crystal_score,
                                                        //  configuration.crystal_tolerance,
                                                         configuration.autostop,
                                                         converged_ligands_p,
                                                         best_so_far_p,
                                                         for_how_long_best_p,
                                                         configuration.tolerance_window,
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
      score_stage.get()->prepare(batch);

    };
    void operator()() {
      auto& chromosomes_b              = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      chromosome* current_population_p = chromosomes_b.dev_pointer();
      chromosome* next_population_p    = next_population.dev_pointer();

      assert(kernel && "Kernel method not yet prepared");
      kernel->set_population_buffers(current_population_p, next_population_p);
      geom_trans.set_chromosomes_buffer(current_population_p);
      kernel->initialize();

      printf("Running GA...\n");
      for (int generation = 0; generation < num_generations; ++generation) {
        geom_trans();
        (*score_stage)();
        (*kernel)();

        // Avoid full device-to-device copy by ping-ponging population buffers.
        if (generation + 1 < num_generations) {
          std::swap(current_population_p, next_population_p);
          kernel->set_population_buffers(current_population_p, next_population_p);
          geom_trans.set_chromosomes_buffer(current_population_p);
        }
      }
      kernel->set_population_buffers(current_population_p, next_population_p);
      kernel->finalize();
    };

    static std::size_t get_shared_ligand_mem(const int max_atoms, const knobs conf) {
      const int chromosomes_per_ligand = std::max(1, static_cast<int>(conf.population_number));
      std::size_t mem{0};
      mem += sizeof(int);                                              // num_atoms
      mem += sizeof(int);                                              // num_rotamers
      mem += sizeof(int);                                              // converged_ligands
      mem += sizeof(chromosome) * chromosomes_per_ligand;              // chromosomes
      mem += sizeof(fp_type) * chromosomes_per_ligand;                 // scores
      mem += 3 * sizeof(fp_type) * max_atoms;                          // coords
      mem += 3 * sizeof(fp_type) * max_atoms * chromosomes_per_ligand; // coord scratch
      return mem;
    }

    static std::size_t get_private_ligand_mem(const int max_atoms, const knobs conf) {
      std::size_t mem{0};
      mem += sizeof(chromosome) * std::max(1, static_cast<int>(conf.population_number)); // next population
      mem += sizeof(chromosome);                                                         // best chromosomes
      mem += sizeof(fp_type);                                                            // best scores
      mem += sizeof(int);                                                                // converged ligands
      mem += scoring_t<queue_t>::get_private_ligand_mem(max_atoms, conf);
      mem += geometric<queue_t>::get_private_ligand_mem(max_atoms, conf);
      return mem;
    }

    static int get_ligand_mem(const int max_atoms, const knobs conf) {
      return static_cast<int>(get_shared_ligand_mem(max_atoms, conf) +
                              get_private_ligand_mem(max_atoms, conf));
    }

    static batch_multiple get_batch_size(const int atoms,
                                         std::shared_ptr<queue_t> q,
                                         const knobs& conf,
                                         const size_t max_bucket_size) {
      (void) max_bucket_size;
      const auto score_bucket_info =
          normalize_batch_multiple(scoring_t<queue_t>::get_batch_size(atoms, q, conf, max_bucket_size));
      const auto geom_bucket_info =
          normalize_batch_multiple(geometric<queue_t>::get_batch_size(atoms, q, conf, max_bucket_size));
      const int score_total = score_bucket_info.total_multiple();
      const int geom_total  = geom_bucket_info.total_multiple();

      batch_multiple selected_info{};
      const char* combine_policy = "MIN";
  #ifdef MUDOCK_GENETIC_BUCKET_COMBINE_SCORE_ONLY
      selected_info  = score_bucket_info;
      combine_policy = "SCORE_ONLY";
  #elif defined(MUDOCK_GENETIC_BUCKET_COMBINE_GEOM_ONLY)
      selected_info  = geom_bucket_info;
      combine_policy = "GEOM_ONLY";
  #elif defined(MUDOCK_GENETIC_BUCKET_COMBINE_LCM)
      {
        const long long lcm_total =
            std::lcm(static_cast<long long>(score_total), static_cast<long long>(geom_total));
        if (lcm_total <= 0 || lcm_total > static_cast<long long>(std::numeric_limits<int>::max())) {
          throw std::runtime_error("GENETIC stage LCM combine overflowed int range");
        }
        // LCM is a pure combined multiplicity; represent it as total x 1.
        selected_info = batch_multiple{static_cast<int>(lcm_total), 1};
      }
      combine_policy = "LCM";
  #else
      if (score_total <= geom_total) {
        selected_info = score_bucket_info;
      } else {
        selected_info = geom_bucket_info;
      }
  #endif
      selected_info = normalize_batch_multiple(selected_info);
      mudock::stage_bucket_trace("GENETIC stage combine for ",
                                 atoms,
                                 " atoms: score_multiple=",
                                 score_total,
                                 " (",
                                 score_bucket_info.active_blocks_per_sm,
                                 "x",
                                 score_bucket_info.num_sms,
                                 "), geom_multiple=",
                                 geom_total,
                                 " (",
                                 geom_bucket_info.active_blocks_per_sm,
                                 "x",
                                 geom_bucket_info.num_sms,
                                 ") policy=",
                                 combine_policy,
                                 " -> selected_plain_multiple=",
                                 selected_info.total_multiple(),
                                 " (",
                                 selected_info.active_blocks_per_sm,
                                 "x",
                                 selected_info.num_sms,
                                 ")");
      return selected_info;
    }
    
  protected:
   
    std::unique_ptr<genetic_kernel<queue_t>> kernel;
    std::shared_ptr<scoring_t<queue_t>> score_stage;
    geometric<queue_t> geom_trans;

    int batch_ligands;
    int num_generations;
    int population_number;
    buffer_vector<chromosome, queue_t> next_population;
    buffer_vector<chromosome, queue_t> best_chromosomes;
    buffer_vector<fp_type, queue_t> best_scores;
    // buffer_vector<int, queue_t> converged_ligands;

    void teardown_impl(batch<static_molecule>& batch) {
      assert(batch.num_ligands == batch_ligands && "Genetic algorithm received different batch for teardown");
  
      // TODO check if the copy can be changed with a swap
      // auto& population_b = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      // population_b.copy_device2device(best_chromosomes);
      // geom_trans();
      // score_stage();
      // geom_trans.teardown(batch);
      // score_stage.teardown(batch);
      auto &converged_ligands_b = (*this->scratch).template get<buffer_data_type::CONVERGED_LIGANDS>();
      best_scores.copy_device2host();
      converged_ligands_b.copy_device2host();
      (*this->scratch).get_queue()->synchronize();
  
      for (int index{0}; index < batch_ligands; ++index) {
        auto& ligand = *batch.molecules[index];
        ligand.properties.assign(property_type::SCORE, std::to_string(best_scores()[index]));

        const int convergence_generation = converged_ligands_b()[index];
        int past_generations = 0;
        if (convergence_generation != 0){
          past_generations = convergence_generation;
        } else {
          past_generations = num_generations;
        }
        const int num_evaluations = past_generations * population_number;
        ligand.properties.assign(property_type::GEN, std::to_string(past_generations));
        ligand.properties.assign(property_type::NUM_EVALS, std::to_string(num_evaluations));
      }
    }

  }; // namespace mudock
#endif
} // namespace mudock
