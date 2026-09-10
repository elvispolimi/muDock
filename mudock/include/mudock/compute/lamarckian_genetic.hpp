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
      genetic<queue_t, scoring_t>::prepare(batch);
      const knobs& configuration = (*this->scratch).configuration;
      ls_every                   = static_cast<int>(configuration.ls_every);
      ls_last_gen                = static_cast<int>(configuration.ls_last_gen);
      local_search_rate          = static_cast<int>(configuration.lsrate);
      local_search_iterations    = static_cast<int>(configuration.lsit);
      
      // Lorenzo: Ligand properties for experiments
      for (int index{0}; index < this->batch_ligands; ++index) {
        auto& ligand = *batch.molecules[index];
        ligand.properties.assign(property_type::RATE, std::to_string(local_search_rate));
        ligand.properties.assign(property_type::ITER, std::to_string(local_search_iterations));
      }
      
      local_search_stage.prepare(batch);  

    };

    void operator()() override {
      auto& chromosomes_b = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      chromosome* current_population_p = chromosomes_b.dev_pointer();
      chromosome* next_population_p    = this->next_population.dev_pointer();
      
      assert(this->kernel && "lamarckian_kernel method not yet prepared");
      this->kernel->set_population_buffers(current_population_p, next_population_p);
      this->geom_trans.set_chromosomes_buffer(current_population_p);
      this->kernel->initialize();

      printf("Running LGA...\n");
      for (int generation = 1; generation <= this->num_generations; ++generation) {
        this->geom_trans();
        (*this->score_stage)();
        (*this->kernel)();

        // Avoid full device-to-device copy by ping-ponging population buffers.
        if (generation < this->num_generations) {
          std::swap(current_population_p, next_population_p);
          this->kernel->set_population_buffers(current_population_p, next_population_p);
          this->geom_trans.set_chromosomes_buffer(current_population_p);
        }

        // Limit LS runs: balance speed and results
        if (generation % ls_every == 0 || generation > (this->num_generations - ls_last_gen)) {
          local_search_stage();
        }

      }
      this->kernel->set_population_buffers(current_population_p, next_population_p);
      this->kernel->finalize();
    }

    static std::size_t get_shared_ligand_mem(const int max_atoms, const knobs conf) {
      std::size_t mem{0};
      mem += genetic<queue_t, scoring_t>::get_shared_ligand_mem(max_atoms, conf);
      return mem;
    }

    static std::size_t get_private_ligand_mem(const int max_atoms, const knobs conf) {
      std::size_t mem{0};
      mem += genetic<queue_t, scoring_t>::get_private_ligand_mem(max_atoms, conf);
      mem += local_search_t<queue_t, scoring_t>::get_private_ligand_mem(max_atoms, conf);
      return mem;
    }

    static int get_ligand_mem(const int max_atoms, const knobs conf) {
      return static_cast<int>(get_shared_ligand_mem(max_atoms, conf) +
                              get_private_ligand_mem(max_atoms, conf));
    }

    void teardown_impl(batch<static_molecule>& batch) override {
      assert(batch.num_ligands == batch_ligands && "Genetic algorithm received different batch for teardown");

      auto &converged_ligands_b = (*this->scratch).template get<buffer_data_type::CONVERGED_LIGANDS>();
      this->best_scores.copy_device2host();
      converged_ligands_b.copy_device2host();
      (*this->scratch).get_queue()->synchronize();
  
      for (int index{0}; index < this->batch_ligands; ++index) {
        auto& ligand = *batch.molecules[index];
        ligand.properties.assign(property_type::SCORE, std::to_string(this->best_scores()[index]));

        const int convergence_generation = converged_ligands_b()[index];
        int past_generations = 0;
        if (convergence_generation != 0){
          past_generations = convergence_generation;
        } else {
          past_generations = this->num_generations;
        }
        // TODO L check if this estimate is correct. WARNING: this depends on the local search implementation. 
        // This is for adadelta for example (not counting the effect of ls_last_gen and ls_every)
        // +1 comes from the scores, the remaining from the gradients
        // TODO L this is not correct: if autostop is on, it should count the actual number of generations at convergence.
        // It would be better to have a counter at each evaluation to be sure (pay attention to race conditions)
        const int num_local_search_individuals = this->population_number * local_search_rate / 100;
        const int num_evaluations = past_generations * this->population_number + past_generations * num_local_search_individuals * local_search_iterations; // GA contribution + LS contribution
        ligand.properties.assign(property_type::GEN, std::to_string(past_generations));
        ligand.properties.assign(property_type::NUM_EVALS, std::to_string(num_evaluations));
      }
    }
  private:
    local_search_t<queue_t, scoring_t> local_search_stage;
    int ls_every;
    int ls_last_gen;
    int local_search_rate;
    int local_search_iterations;
  }; // namespace mudock
#endif
} // namespace mudock
