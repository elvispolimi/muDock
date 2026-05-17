#pragma once

#include <cstddef>
#include <cstring>
#include <functional>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/compute/adadelta_kernel.hpp>
#include <mudock/compute/local_search.hpp>
#ifndef __CUDACC__
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  #define ADADELTA_RHO 0.85f
  #define ADADELTA_EPSILON 1e-6f
  #define ADADELTA_CONVERGENCE_THRESHOLD 1e-6f
  #define ADADELTA_CONVERGENCE_PATIENCE 5

  #ifndef __CUDACC__
  // TODO check that the object type and the kernel impl are the same
  template<typename queue_type, template<typename> typename scoring_t>
  struct adadelta: public local_search<queue_type, scoring_t> {
    static constexpr const char stage_name[] = "ADADELTA";
    
    adadelta(std::shared_ptr<scratchpad<queue_type>> _scratch,
             std::shared_ptr<scoring_t<queue_type>> _score) 
             : local_search<queue_type, scoring_t>(_scratch, _score) {}

    void prepare(batch<static_molecule> &batch) {
      // TODO L i don't like initializing iterations here, not scalable. Better move it to local_search
      this->iterations  = (*this->scratch).configuration.ls_iterations;
      // Allocate gradient buffer for AdaDelta (one gradient per individual per ligand)
      batch_ligands = batch.num_ligands;
      const int individuals_per_ligand = std::max(1, static_cast<int>((*this->scratch).configuration.population_number));
      auto &gradient_b = (*this->scratch).template get<buffer_data_type::GRADIENTS>();
      const size_t gradient_count = static_cast<size_t>(batch_ligands) * static_cast<size_t>(individuals_per_ligand);
      if (!gradient_b.is_valid() || gradient_b.num_elements() != gradient_count) {
        gradient_b.alloc(gradient_count);
        gradient_b.set_valid();
      }
      gradient *gradients_b = gradient_b.dev_pointer();

      auto &chromosomes_b = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      chromosomes_b.alloc(static_cast<size_t>(batch_ligands) * individuals_per_ligand);
      chromosomes_b.set_valid();
      chromosome *population_b = chromosomes_b.dev_pointer();

      // Allocate AdaDelta state buffers (E[g^2] and E[delta^2])
      auto &adadelta_e_g2_b = (*this->scratch).template get<buffer_data_type::ADADELTA_E_G2>();
      auto &adadelta_e_dw2_b = (*this->scratch).template get<buffer_data_type::ADADELTA_E_DW2>();
      auto &stall_counter_b = (*this->scratch).template get<buffer_data_type::STALL_COUNTER>();
      auto &inactive_b = (*this->scratch).template get<buffer_data_type::INACTIVE>();

      if (!adadelta_e_g2_b.is_valid() || adadelta_e_g2_b.num_elements() != gradient_count) {
        adadelta_e_g2_b.alloc(gradient_count);
        adadelta_e_dw2_b.alloc(gradient_count);
        stall_counter_b.alloc(gradient_count);
        inactive_b.alloc(gradient_count);

        adadelta_e_g2_b.set_valid();
        adadelta_e_dw2_b.set_valid();
        stall_counter_b.set_valid();
        inactive_b.set_valid();

      }
      chromosome *adadelta_e_g2 = adadelta_e_g2_b.dev_pointer();
      chromosome *adadelta_e_dw2 = adadelta_e_dw2_b.dev_pointer();
      int *stall_counter = stall_counter_b.dev_pointer();
      int *inactive = inactive_b.dev_pointer();

      auto q = (*this->scratch).get_queue();

      ls_ad_kernel = std::make_unique<adadelta_kernel<queue_type>>(individuals_per_ligand,
                                                                  batch_ligands,
                                                                  this->score_stage,
                                                                  gradients_b,
                                                                  population_b,
                                                                  adadelta_e_g2,
                                                                  adadelta_e_dw2,
                                                                  stall_counter,
                                                                  inactive,
                                                                  q);

      // Initialize the scoring kernel buffers
      this->score_stage->prepare(batch);
      
    }

    void operator()() {
      assert(
          (((*this->scratch).template get<buffer_data_type::GRADIENTS>().num_elements() % batch_ligands) == 0) &&
          "Number of gradients is not a multiple of ligands in the batch");
      assert(ls_ad_kernel && "Adadelta local search kernel method not yet prepared");

      // TODO L try to move the reset of inactives here which is more elegant, for now it is in adadelta cpp


      for (int i = 0; i < this->iterations; ++i) {
        if (coordinate_update) {
          coordinate_update();
        }
        ls_ad_kernel->compute_gradients();
        ls_ad_kernel->apply_adadelta(i);
      }

    }

    void set_coordinate_update(std::function<void()> update) { coordinate_update = std::move(update); }

    static int get_ligand_mem(const int max_atoms, const knobs conf) {
      int mem{0};
      const int individuals_per_ligand = std::max(1, static_cast<int>(conf.population_number));

      // One gradient per individual per ligand
      mem += sizeof(gradient) * individuals_per_ligand;
      // AdaDelta state buffers for each individual
      mem += sizeof(chromosome) * individuals_per_ligand; // E[g^2]
      mem += sizeof(chromosome) * individuals_per_ligand; // E[delta_w^2]
      mem += sizeof(int) * individuals_per_ligand;        // stall counter
      mem += sizeof(int) * individuals_per_ligand;        // inactive flag
      return mem;
    }

    static batch_multiple get_batch_size(const int atoms,
                                         std::shared_ptr<queue_type> q,
                                         const knobs &conf,
                                         const size_t max_bucket_size) {
      (void) conf;
      const auto plain_multiple_info =
          normalize_batch_multiple(get_adt_score_batch_multiple<queue_type>(atoms, q));
      mudock::stage_bucket_trace("ADADELTA stage plain multiple for ",
                                 atoms,
                                 " atoms -> total=",
                                 plain_multiple_info.total_multiple(),
                                 " (active_blocks_per_sm=",
                                 plain_multiple_info.active_blocks_per_sm,
                                 ", num_sms=",
                                 plain_multiple_info.num_sms,
                                 ")",
                                 " (max_bucket_size hint=",
                                 max_bucket_size,
                                 ")");
      return plain_multiple_info;
    }

  private:
    int batch_ligands;

    std::unique_ptr<adadelta_kernel<queue_type>> ls_ad_kernel;
    std::function<void()> coordinate_update;

    void teardown_impl(batch<static_molecule> &batch) override {
      // TODO L: Implement teardown
    };
  };
  #endif
} // namespace mudock
