#pragma once

#include "mudock/knobs.hpp"
#include "genetic.hpp"

#include <concepts>
#include <memory>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_protein.hpp>
#if !defined(__CUDACC__) && !defined(__HIPCC__)
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/docking.hpp>
  #include <mudock/compute/geometric_transform.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/local_search.hpp>
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
                              const fp_type mutation_prob_,
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
    inline void set_population_buffers(chromosome* population_, chromosome* next_population_) { genetic_kernel<queue_type>::set_population_buffers(population_, next_population_); }

  private:
  };

#if !defined(__CUDACC__) && !defined(__HIPCC__)
  template<typename queue_t, template<typename> typename scoring_t, template<typename, template<typename> typename> typename local_search_t>
    requires std::derived_from<queue_t, queue> 
            && std::derived_from<scoring_t<queue_t>, scoring<queue_t>>
            && std::derived_from<local_search_t<queue_t, scoring_t>, local_search<queue_t, scoring_t>>
  struct lamarckian_genetic : public genetic<queue_t, scoring_t> {
    static constexpr const char stage_name[] = "LAMARCKIAN_GENETIC";
    // Inherit constructor and most logic from genetic
    lamarckian_genetic(std::shared_ptr<scratchpad<queue_t>> _scratch,
                      dynamic_molecule& _protein,
                      std::shared_ptr<scoring_t<queue_t>> _scoring,
                      local_search_t<queue_t, scoring_t> _local_search)
        : genetic<queue_t, scoring_t>(_scratch, _protein, _scoring),
          local_search_stage(std::move(_local_search)){};

    void prepare(batch<static_molecule>& batch) override {
      const knobs& configuration  = (*this->scratch).configuration;
      this->batch_ligands         = batch.num_ligands;
      this->num_generations       = static_cast<int>(configuration.num_generations);
      const int population_number = static_cast<int>(configuration.population_number);
      auto q                      = (*this->scratch).get_queue();

      auto& num_rotamers_b = (*this->scratch).template get<buffer_data_type::NUM_ROTAMERS>();
      auto& chromosomes_b  = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      auto& scores_b       = (*this->scratch).template get<buffer_data_type::SCORES>();

      num_rotamers_b.alloc(this->batch_ligands);
      chromosomes_b.alloc(population_number * this->batch_ligands);
      this->next_population.alloc(population_number * this->batch_ligands);
      scores_b.alloc(population_number * this->batch_ligands);
      this->best_scores.alloc(this->batch_ligands);
      this->best_chromosomes.alloc(this->batch_ligands);

      load_num_rotamers<queue_t>(batch, this->scratch);

      const auto seed =
          configuration.seed.has_value()
              ? configuration.seed.value()
              : static_cast<size_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());

      int* __restrict__ num_rotamers_p            = num_rotamers_b.dev_pointer();
      fp_type* __restrict__ scores_p              = scores_b.dev_pointer();
      fp_type* __restrict__ best_scores_p         = this->best_scores.dev_pointer();
      chromosome* __restrict__ best_chromosomes_p = this->best_chromosomes.dev_pointer();

      this->lamarckian_kernel = std::make_unique<lamarckian_genetic_kernel<queue_t>>(this->batch_ligands,
                                                                                    population_number,
                                                                                    configuration.num_generations,
                                                                                    configuration.tournament_length,
                                                                                    configuration.mutation_prob,
                                                                                    seed,
                                                                                    chromosomes_b.dev_pointer(),
                                                                                    this->next_population.dev_pointer(),
                                                                                    num_rotamers_p,
                                                                                    scores_p,
                                                                                    best_scores_p,
                                                                                    best_chromosomes_p,
                                                                                    q);

      this->geom_trans.prepare(batch);
      local_search_stage.prepare(batch);  
      (this->score_stage).get()->prepare(batch);
      // TODO L at the moment this prepare() call order must be kept
      // (ls and then score) in order to initialize correctly active population also in score. Try to make it independent
      // maybe moving it here in LGA
    };

    void operator()() override {
      auto& chromosomes_b = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      chromosome* current_population_p = chromosomes_b.dev_pointer();
      chromosome* next_population_p    = this->next_population.dev_pointer();
      
      assert(this->lamarckian_kernel && "lamarckian_kernel method not yet prepared");
      this->lamarckian_kernel->set_population_buffers(current_population_p, next_population_p);
      this->geom_trans.set_chromosomes_buffer(current_population_p);
      this->lamarckian_kernel->initialize();

      printf("Running LGA...\n");
      for (int generation = 0; generation < this->num_generations; ++generation) {
        this->geom_trans();
        (*this->score_stage)();
        local_search_stage();
        (*this->lamarckian_kernel)();

        // Avoid full device-to-device copy by ping-ponging population buffers.
        if (generation + 1 < this->num_generations) {
          std::swap(current_population_p, next_population_p);
          this->lamarckian_kernel->set_population_buffers(current_population_p, next_population_p);
          this->geom_trans.set_chromosomes_buffer(current_population_p);
        }
      }
      this->lamarckian_kernel->set_population_buffers(current_population_p, next_population_p);
      this->lamarckian_kernel->finalize();
    }

    // TODO L Important: this was copied from genetic.hpp code, but not sure if it must be adapted 
    static std::size_t get_shared_ligand_mem(const int max_atoms, const knobs conf) {
      const int chromosomes_per_ligand = std::max(1, static_cast<int>(conf.population_number));
      std::size_t mem{0};
      mem += sizeof(int);                                              // num_atoms
      mem += sizeof(int);                                              // num_rotamers
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
      mem += scoring_t<queue_t>::get_private_ligand_mem(max_atoms, conf);
      mem += geometric<queue_t>::get_private_ligand_mem(max_atoms, conf);
      mem += local_search_t<queue_t, scoring_t>::get_private_ligand_mem(max_atoms, conf);
      return mem;
    }

    static int get_ligand_mem(const int max_atoms, const knobs conf) {
      return static_cast<int>(get_shared_ligand_mem(max_atoms, conf) +
                              get_private_ligand_mem(max_atoms, conf));
    }

    // TODO L Important: this was copied from genetic.hpp code, but not sure if it must be adapted 
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
          throw std::runtime_error("LAMARCKIAN GENETIC stage LCM combine overflowed int range");
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
      mudock::stage_bucket_trace("LAMARCKIAN GENETIC stage combine for ",
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
    local_search_t<queue_t, scoring_t> local_search_stage;
    std::unique_ptr<lamarckian_genetic_kernel<queue_t>> lamarckian_kernel;
    
    void teardown_impl(batch<static_molecule>& batch) {
      assert(batch.num_ligands == this->batch_ligands && "Lamarckian-Genetic algorithm received different batch for teardown");

      this->best_scores.copy_device2host();
      (*this->scratch).get_queue()->synchronize();

      for (int index{0}; index < this->batch_ligands; ++index) {
        auto& ligand = *batch.molecules[index];
        ligand.properties.assign(property_type::SCORE, std::to_string(this->best_scores()[index]));
      }
    }
  }; // namespace mudock
#endif
} // namespace mudock
