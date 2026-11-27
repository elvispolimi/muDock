#pragma once

#include "mudock/knobs.hpp"

#include <algorithm>
#include <chrono>
#include <concepts>
#include <cstddef>
#include <memory>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <utility>
#include <mudock/compute/buffer_utils.hpp>
#include <mudock/compute/docking.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/compute/ligand_placement.hpp>
#include <mudock/compute/scoring.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/compute/genetic_args.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/compute/geom_transform_args.hpp>
#include <mudock/compute/algorithm.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  template<typename queue_t>
  struct genetic_backend {
    template<scoring_function scoring>
    static void run(std::shared_ptr<queue_t> queue,
                    const genetic_args& genetic,
                    geometric<queue_t>& geometry_stage,
                    typename scoring_backend<queue_t, scoring>::stage_type& scoring_stage_instance);
  };

  template<typename queue_t, scoring_function scoring>
    requires std::derived_from<queue_t, queue>
  struct genetic: public docking<queue_t> {
    using scoring_stage = typename scoring_backend<queue_t, scoring>::stage_type;
    static constexpr const char stage_name[] = "GENETIC";

    genetic(std::shared_ptr<scratchpad<queue_t>> _scratch,
            dynamic_molecule& protein,
            scoring_stage _scoring,
            ligand_placement _placement = {})
        : docking<queue_t>(_scratch),
          score_stage(std::move(_scoring)),
          geom_trans(_scratch, protein, _placement),
          next_population(_scratch->get_queue()),
          best_chromosomes(_scratch->get_queue()),
          best_scores(_scratch->get_queue()) {};
    void prepare(batch<static_molecule>& batch) {
      trace_id = batch.trace_id;
      trace_event("batch_received");
      const knobs& configuration  = (*this->scratch).configuration;
      batch_ligands               = batch.num_ligands;
      const int population_number = static_cast<int>(configuration.population_number);

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

      args = {batch_ligands,
              population_number,
              static_cast<int>(configuration.num_generations),
              seed,
              static_cast<int>(configuration.tournament_length),
              configuration.mutation_prob,
              configuration.max_translation_half_width.value_or(fp_type{0}),
              chromosomes_b.dev_pointer(),
              next_population.dev_pointer(),
              num_rotamers_p,
              scores_p,
              best_scores_p,
              best_chromosomes_p,
              {}};

      geom_trans.prepare(batch);
      score_stage.prepare(batch);
      args.geometry = geom_trans.kernel_args();
      trace_event("input_copies_submitted");
    };
    void operator()() override {
      trace_event("compute_enqueue_begin");
      genetic_backend<queue_t>::template run<scoring>(
          this->scratch->get_queue(), args, geom_trans, score_stage);
      trace_event("compute_enqueue_end");
    }

    [[nodiscard]] genetic_args kernel_args() const {
      return args;
    }

    static std::size_t get_shared_ligand_mem(const int max_atoms, const knobs conf) {
      const int chromosomes_per_ligand = std::max(1, static_cast<int>(conf.population_number));
      std::size_t mem{0};
      mem += sizeof(int);                                              // num_atoms
      mem += sizeof(int);                                              // num_rotamers
      mem += sizeof(chromosome) * chromosomes_per_ligand;              // chromosomes
      mem += sizeof(fp_type) * chromosomes_per_ligand;                 // scores
      mem += 3 * sizeof(fp_type) * max_atoms;                          // coords
      return mem;
    }

    static std::size_t get_private_ligand_mem(const int max_atoms, const knobs conf) {
      std::size_t mem{0};
      mem += sizeof(chromosome) * std::max(1, static_cast<int>(conf.population_number)); // next population
      mem += sizeof(chromosome);                                                         // best chromosomes
      mem += sizeof(fp_type);                                                            // best scores
      mem += scoring_stage::get_private_ligand_mem(max_atoms, conf);
      mem += geometric<queue_t>::get_private_ligand_mem(max_atoms, conf);
      return mem;
    }

    static std::size_t get_ligand_mem(const int max_atoms, const knobs conf) {
      return get_shared_ligand_mem(max_atoms, conf) + get_private_ligand_mem(max_atoms, conf);
    }

  private:
    std::size_t trace_id = 0;
    scoring_stage score_stage;
    geometric<queue_t> geom_trans;
    genetic_args args{};

    int batch_ligands;
    buffer_vector<chromosome, queue_t> next_population;
    buffer_vector<chromosome, queue_t> best_chromosomes;
    buffer_vector<fp_type, queue_t> best_scores;

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

      trace_event("copyback_begin");
      best_scores.copy_device2host();
      best_chromosomes.copy_device2host();
      (*this->scratch).get_queue()->synchronize();
      trace_event("compute_and_copyback_complete");

      for (int index{0}; index < batch_ligands; ++index) {
        auto& ligand = *batch.molecules[index];
        geom_ligand geom{ligand};
        const int num_rotamers = static_cast<int>(ligand.num_rotamers());
        apply<cpu_vectorization::AUTO>(ligand.x(),
                                       ligand.y(),
                                       ligand.z(),
                                       best_chromosomes()[index],
                                       ligand.num_atoms(),
                                       num_rotamers,
                                       geom.fragments_masks(),
                                       geom.fragmets_starts(),
                                       geom.fragments_stops());
        ligand.properties.assign(property_type::SCORE, std::to_string(best_scores()[index]));
      }
    }
  }; // namespace mudock
} // namespace mudock
