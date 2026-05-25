#pragma once

#include <cstddef>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <optional>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/compute/adadelta_kernel.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/compute/local_search.hpp>
#ifndef __CUDACC__
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  #define RHO 0.8f
  #define EPSILON 1e-2f                // TODO L use variables instead of numbers for coordinate and angle step
  #define CONVERGENCE_THRESHOLD_COORD (static_cast<fp_type>(0.2 * 90 * 0.01)) // 90*0.2 = 18 is the range of movement for translations
  #define CONVERGENCE_THRESHOLD_ANGLE (static_cast<fp_type>(4 * 90 * 0.01))      // 90*4 = 360 is the range of movement for angles
  #define MAX_STEP 1e-1f
  #define MAX_STEP_POS  0.1f
  #define MAX_STEP_ROT  0.03f
  #define MAX_STEP_TORS 0.05f

  #ifndef __CUDACC__
  // TODO check that the object type and the kernel impl are the same
  template<typename queue_type, template<typename> typename scoring_t>
  struct adadelta: public local_search<queue_type, scoring_t> {
    static constexpr const char stage_name[] = "ADADELTA";
    
    adadelta(std::shared_ptr<scratchpad<queue_type>> _scratch,
             std::shared_ptr<scoring_t<queue_type>> _score) 
             : local_search<queue_type, scoring_t>(_scratch, _score),
               geom_trans(_scratch, _score->get_protein()) {}

    void prepare(batch<static_molecule> &batch) {
      // TODO L check if this makes a copy or a reference
      assert(batch.num_ligands == 1 && "AdaDelta dump_pose currently expects a single ligand in the batch");
      this->ligand_template = *batch.molecules[0];

      batch_ligands = batch.num_ligands;
      const int individuals_per_ligand = std::max(1, static_cast<int>((*this->scratch).configuration.population_number));
      auto &gradient_b = (*this->scratch).template get<buffer_data_type::GRADIENTS>();
      const size_t gradient_count = static_cast<size_t>(batch_ligands) * static_cast<size_t>(individuals_per_ligand);
      if (!gradient_b.is_valid() || gradient_b.num_elements() != gradient_count) {
        gradient_b.alloc(gradient_count);
        gradient_b.set_valid();
      }
      gradient *gradients_b = gradient_b.dev_pointer();

      auto& num_rotamers_b = (*this->scratch).template get<buffer_data_type::NUM_ROTAMERS>();
      num_rotamers_b.alloc(batch_ligands);
      load_num_rotamers<queue_type>(batch, this->scratch);

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

      int* __restrict__ num_rotamers_p            = num_rotamers_b.dev_pointer();
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
                                                                  num_rotamers_p,
                                                                  adadelta_e_g2,
                                                                  adadelta_e_dw2,
                                                                  stall_counter,
                                                                  inactive,
                                                                  q);

      // Initialize the scoring kernel buffers
      this->score_stage->prepare(batch);
      geom_trans.prepare(batch);
    }

    void operator()() {
      assert(
          (((*this->scratch).template get<buffer_data_type::GRADIENTS>().num_elements() % batch_ligands) == 0) &&
          "Number of gradients is not a multiple of ligands in the batch");
      assert(ls_ad_kernel && "Adadelta local search kernel method not yet prepared");

      auto &scores_b = (*this->scratch).template get<buffer_data_type::SCORES>();

      // TODO L try to move the reset of inactives here which is more elegant, for now it is in adadelta cpp

      const bool only_local_search =
          ((*this->scratch).configuration.population_number == 1) &&
          ((*this->scratch).configuration.num_generations == 1);

      for (std::size_t i = 0; i < this->iterations; ++i) {
        geom_trans();
        
        if (only_local_search) {
          (*this->score_stage)();

          // copy scores back to host and print best score (first element)
          scores_b.copy_device2host();
          (*this->scratch).get_queue()->synchronize();
          if(i % (this->iterations/10) == 0){ // print eveery 10% of the process
            this->dump_pose(int(i));
            printf("Iter: %ld, Score: %f\n", i, scores_b()[0]);
          }
        }

        ls_ad_kernel->compute_gradients();
        ls_ad_kernel->apply_adadelta(static_cast<int>(i), static_cast<int>(this->convergence_patience), this->use_early_stopping);
      }

      if (only_local_search) {
        geom_trans();
        (*this->score_stage)();
      }
    }

    static int get_ligand_mem(const int max_atoms, const knobs conf) {
      std::size_t mem{0};
      const int individuals_per_ligand = std::max(1, static_cast<int>(conf.population_number));

      mem += scoring_t<queue_type>::get_ligand_mem(max_atoms, conf);

      // One gradient per individual per ligand
      mem += sizeof(gradient) * individuals_per_ligand;
      // AdaDelta state buffers for each individual
      mem += sizeof(chromosome) * individuals_per_ligand; // E[g^2]
      mem += sizeof(chromosome) * individuals_per_ligand; // E[delta_w^2]
      mem += sizeof(int) * individuals_per_ligand;        // stall counter
      mem += sizeof(int) * individuals_per_ligand;        // inactive flag
      return static_cast<int>(mem);
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
    geometric<queue_type> geom_trans;
    std::function<void()> coordinate_update;

    // TODO L: Implement teardown
    void teardown_impl(batch<static_molecule> &batch) override {
      auto &scores_b              = (*this->scratch).template get<buffer_data_type::SCORES>();
      const int scores_per_ligand = static_cast<int>(scores_b.num_elements() / batch_ligands);
      scores_b.copy_device2host();
      (*this->scratch).get_queue()->synchronize();
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand          = *batch.molecules[ligand_index];
        const int score_index = ligand_index * scores_per_ligand;
        ligand.properties.assign(property_type::SCORE, std::to_string(scores_b()[score_index]));
      }
    };

  };
  #endif
} // namespace mudock
