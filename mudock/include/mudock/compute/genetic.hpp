#pragma once

#include "mudock/knobs.hpp"

#include <algorithm>
#include <chrono>
#include <concepts>
#include <cstddef>
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
  #include <mudock/compute/ligand_placement.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/data_formats/data_store.hpp>
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
                   const int tournament_length_,
                   const fp_type mutation_prob_,
                   const fp_type max_translation_half_width_,
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
          max_translation_half_width(max_translation_half_width_),
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
    int tournament_length;
    fp_type mutation_prob;
    fp_type max_translation_half_width;
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
  template<typename queue_t,
           template<typename, typename> typename scoring_t,
           typename layout_set = genetic_data_layouts>
    requires std::derived_from<queue_t, queue> &&
             std::derived_from<scoring_t<queue_t, layout_set>, scoring<queue_t, layout_set>>
  struct genetic: public docking<queue_t, layout_set> {
    static constexpr const char stage_name[] = "GENETIC";

    genetic(std::shared_ptr<data_store<queue_t, layout_set>> _data,
            dynamic_molecule& protein,
        scoring_t<queue_t, layout_set> _scoring,
        ligand_placement _placement = {})
      : docking<queue_t, layout_set>(_data),
          score_stage(std::move(_scoring)),
          geom_trans(_data, protein, _placement) {};
    void prepare(batch<static_molecule>& batch) {
      trace_id = batch.trace_id;
      trace_event("batch_received");
      const knobs& configuration  = (*this->data).configuration;
      batch_ligands               = batch.num_ligands;
      num_generations             = static_cast<int>(configuration.num_generations);
      const int population_number = static_cast<int>(configuration.population_number);
      auto q                      = (*this->data).get_queue();

      const int batch_atoms = batch.batch_max_atoms;
      const int total_atoms = batch_ligands * batch_atoms;
      this->data->alloc(total_atoms,
                           total_atoms * population_number,
                           batch_ligands,
                           population_number * batch_ligands,
                           population_number * batch_ligands,
                           batch_ligands);

      load_num_rotamers<queue_t>(batch, this->data);

      const auto seed =
          configuration.seed.has_value()
              ? configuration.seed.value()
              : static_cast<size_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());

      auto batch_view = this->data->device_view().batch_info();
      auto genetic_view = this->data->device_view().genetic();
      auto scores_view = this->data->device_view().scores();
      auto best_view = this->data->device_view().best();
      int* __restrict__ num_rotamers_p            = batch_view.num_rotamers().data();
      fp_type* __restrict__ scores_p              = scores_view.scores().data();
      fp_type* __restrict__ best_scores_p         = best_view.scores().data();
      chromosome* __restrict__ best_chromosomes_p = best_view.chromosomes().data();

      kernel = std::make_unique<genetic_kernel<queue_t>>(batch_ligands,
                                                         population_number,
                                                         configuration.num_generations,
                                                         configuration.tournament_length,
                                                         configuration.mutation_prob,
                                                         configuration.max_translation_half_width.value_or(fp_type{0}),
                                                         seed,
                                                         genetic_view.chromosomes().data(),
                                                         genetic_view.next_chromosomes().data(),
                                                         num_rotamers_p,
                                                         scores_p,
                                                         best_scores_p,
                                                         best_chromosomes_p,
                                                         q);

      geom_trans.prepare(batch);
      score_stage.prepare(batch);
      trace_event("input_copies_submitted");
    };
    void operator()() {
      auto q                           = (*this->data).get_queue();
      auto genetic_view = this->data->device_view().genetic();
      chromosome* current_population_p = genetic_view.chromosomes().data();
      chromosome* next_population_p    = genetic_view.next_chromosomes().data();

      assert(kernel && "Kernel method not yet prepared");
      kernel->set_population_buffers(current_population_p, next_population_p);
      geom_trans.set_chromosomes_buffer(current_population_p);
      trace_event("compute_enqueue_begin");
      kernel->initialize();
      for (int generation = 0; generation < num_generations; ++generation) {
        geom_trans();
        score_stage();
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
      trace_event("compute_enqueue_end");
    };

    static std::size_t get_shared_ligand_mem(const int max_atoms, const knobs conf) {
      const int chromosomes_per_ligand = std::max(1, static_cast<int>(conf.population_number));
      std::size_t mem{0};
      mem += sizeof(int);                                              // num_atoms
      mem += sizeof(int);                                              // num_rotamers
      mem += sizeof(chromosome) * chromosomes_per_ligand;              // chromosomes
      mem += sizeof(fp_type) * chromosomes_per_ligand;                 // scores
      mem += 3 * sizeof(fp_type) * max_atoms;                          // coords
      mem += 3 * sizeof(fp_type) * max_atoms * chromosomes_per_ligand; // coord data
      return mem;
    }

    static std::size_t get_private_ligand_mem(const int max_atoms, const knobs conf) {
      std::size_t mem{0};
      mem += sizeof(chromosome) * std::max(1, static_cast<int>(conf.population_number)); // next population
      mem += sizeof(chromosome);                                                         // best chromosomes
      mem += sizeof(fp_type);                                                            // best scores
      mem += scoring_t<queue_t, layout_set>::get_private_ligand_mem(max_atoms, conf);
      mem += geometric<queue_t, layout_set>::get_private_ligand_mem(max_atoms, conf);
      return mem;
    }

    static std::size_t get_ligand_mem(const int max_atoms, const knobs conf) {
      return get_shared_ligand_mem(max_atoms, conf) + get_private_ligand_mem(max_atoms, conf);
    }

    static batch_multiple get_batch_size(const int atoms,
                                         std::shared_ptr<queue_t> q,
                                         const knobs& conf,
                                         const size_t max_bucket_size) {
      (void) max_bucket_size;
      const auto score_bucket_info =
          normalize_batch_multiple(scoring_t<queue_t, layout_set>::get_batch_size(atoms, q, conf, max_bucket_size));
      const auto geom_bucket_info =
          normalize_batch_multiple(geometric<queue_t, layout_set>::get_batch_size(atoms, q, conf, max_bucket_size));
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

  private:
    std::size_t trace_id = 0;
    scoring_t<queue_t, layout_set> score_stage;
    geometric<queue_t, layout_set> geom_trans;
    std::unique_ptr<genetic_kernel<queue_t>> kernel;

    int batch_ligands;
    int num_generations;
    void trace_event(const char* event) const {
      if (trace_id == 0)
        return;
      const auto timestamp_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                    std::chrono::steady_clock::now().time_since_epoch())
                                    .count();
      info("Genetic: event=", event, " batch_id=", trace_id, " t_ns=", timestamp_ns);
    }

    void teardown_impl(batch<static_molecule>& batch) {
      assert(batch.num_ligands == batch_ligands && "Genetic algorithm received different batch for teardown");

      // TODO check if the copy can be changed with a swap
      // auto& population_b = (*this->data).template get<buffer_data_type::CHROMOSOMES>();
      // population_b.copy_device2device(best_chromosomes);
      // geom_trans();
      // score_stage();
      // geom_trans.teardown(batch);
      // score_stage.teardown(batch);

      trace_event("copyback_begin");
      this->data->copy_device2host_best();
      (*this->data).get_queue()->synchronize();
      auto best_view = this->data->host_view().best();
      trace_event("compute_and_copyback_complete");

      for (int index{0}; index < batch_ligands; ++index) {
        auto& ligand = *batch.molecules[index];
        geom_ligand geom{ligand};
        const int num_rotamers = static_cast<int>(ligand.num_rotamers());
        apply<cpu_vectorization::AUTO>(ligand.x(),
                                       ligand.y(),
                                       ligand.z(),
                                       best_view.chromosomes()[index],
                                       ligand.num_atoms(),
                                       num_rotamers,
                                       geom.fragments_masks(),
                                       geom.fragmets_starts(),
                                       geom.fragments_stops());
        ligand.properties.assign(property_type::SCORE, std::to_string(best_view.scores()[index]));
      }
    }
  }; // namespace mudock
#endif
} // namespace mudock
