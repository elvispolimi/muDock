#pragma once

#include "mudock/knobs.hpp"

#include <algorithm>
#include <concepts>
#include <limits>
#include <memory>
#include <stdexcept>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <numeric>
#if !defined(__CUDACC__) && !defined(__HIPCC__)
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/docking.hpp>
  #include <mudock/compute/geometric_transform.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/compute/queue.hpp>
#include <mudock/compute/batch_multiple.hpp>
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
          seed(seed_),
          q(q_) {};
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
    size_t seed;
    std::shared_ptr<queue_type> q;
  };

#if !defined(__CUDACC__) && !defined(__HIPCC__)
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
                                                         best_chromosomes_p,
                                                         q);

      geom_trans.prepare(batch);
      score_stage.prepare(batch);
    };
    void operator()() {
      auto& chromosomes_b = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();

      assert(kernel && "Kernel method not yet prepared");
      kernel->initialize();

      for (int generation = 0; generation < num_generations; ++generation) {
        geom_trans();
        score_stage();
        (*kernel)();
        chromosomes_b.copy_device2device(next_population);
      }
      kernel->finalize();
    };

    static int get_ligand_mem(const int max_atoms, const knobs conf) {
      int mem{0};

      mem += sizeof(int);                                 // num_rotamers
      mem += sizeof(chromosome) * conf.population_number; // chromosomes
      mem += sizeof(chromosome) * conf.population_number; //next population
      mem += sizeof(fp_type) * conf.population_number;    //scores
      mem += sizeof(chromosome);                          // best chromosomes
      mem += sizeof(fp_type);                             // best scores

      mem += scoring_t<queue_t>::get_ligand_mem(max_atoms, conf);
      mem += geometric<queue_t>::get_ligand_mem(max_atoms, conf);
      return mem;
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
      selected_info              = score_bucket_info;
      combine_policy             = "SCORE_ONLY";
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
      mudock::info("GENETIC stage combine for ",
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

  private:
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
