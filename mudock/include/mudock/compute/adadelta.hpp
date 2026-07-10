#pragma once

#include <cstddef>
#include <cstring>
#include <random>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/compute/adadelta_kernel.hpp>
#include <mudock/compute/geometric_transform.hpp>
#if !defined(__CUDACC__) && !defined(__HIPCC__)
  #include <mudock/compute/local_search.hpp>
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<typename queue_type>
  batch_multiple get_adadelta_batch_multiple(const int, std::shared_ptr<queue_type>) {
    return {};
  }

#if !defined(__CUDACC__) && !defined(__HIPCC__)
  // TODO check that the object type and the kernel impl are the same
  template<typename queue_type, template<typename> typename scoring_t>
  struct adadelta: public local_search<queue_type, scoring_t> {
    static constexpr const char stage_name[] = "ADADELTA";
    static constexpr fp_type RHO     = 0.8f;
    static constexpr fp_type EPSILON = 1e-2f;
    
    adadelta(std::shared_ptr<scratchpad<queue_type>> _scratch,
             std::shared_ptr<scoring_t<queue_type>> _score) 
             : local_search<queue_type, scoring_t>(_scratch, _score),
               geom_trans(_scratch, _score->get_protein()) {}

    void prepare(batch<static_molecule> &batch) {
      // TODO L check if this makes a copy or a reference
      assert(batch.num_ligands == 1 && "AdaDelta dump_pose currently expects a single ligand in the batch");
      this->ligand_template = *batch.molecules[0];

      batch_ligands = batch.num_ligands;
      batch_atoms   = batch.batch_max_atoms;
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
      chromosomes_b.alloc(static_cast<size_t>(batch_ligands) * static_cast<size_t>(individuals_per_ligand));
      chromosomes_b.set_valid();
      chromosome *population_b = chromosomes_b.dev_pointer();

      // Allocate AdaDelta state buffers (E[g^2] and E[delta^2])
      auto &adadelta_e_g2_b      = (*this->scratch).template get<buffer_data_type::ADADELTA_E_G2>();
      auto &adadelta_e_dw2_b     = (*this->scratch).template get<buffer_data_type::ADADELTA_E_DW2>();
      
      if (!adadelta_e_g2_b.is_valid() || adadelta_e_g2_b.num_elements() != gradient_count) {
        adadelta_e_g2_b.alloc(gradient_count);
        adadelta_e_dw2_b.alloc(gradient_count);
        adadelta_e_g2_b.set_valid();
        adadelta_e_dw2_b.set_valid();
      }
      
      
      // TODO L move active from lsrate initialization from here to local search generic or LGA?
      // Initialize active population
      auto &active_individuals_b = (*this->scratch).template get<buffer_data_type::ACTIVE_INDIVIDUALS>();
      if (!active_individuals_b.is_valid() || active_individuals_b.num_elements() != gradient_count) {
        active_individuals_b.alloc(gradient_count);
        active_individuals_b.set_valid();
      }
      std::vector<int> active_init(gradient_count);
      
      std::mt19937 rng((*this->scratch).configuration.seed.value_or(std::random_device{}()));
      std::bernoulli_distribution dist(static_cast<double>((*this->scratch).configuration.lsrate / fp_type{100}));

      for (size_t i = 0; i < gradient_count; ++i) {
        active_init[i] = dist(rng);
      }
// TODO L IMPORTANT magari devo fare copyhost2device di active individuals??? come in adt_score.hpp
      // Copy to device/managed buffer
      std::memcpy((void *) active_individuals_b(),
                  active_init.data(),
                  gradient_count * sizeof(int));
      active_individuals_b.copy_host2device();

      int* __restrict__ num_rotamers_p = num_rotamers_b.dev_pointer();
      chromosome *adadelta_e_g2        = adadelta_e_g2_b.dev_pointer();
      chromosome *adadelta_e_dw2       = adadelta_e_dw2_b.dev_pointer();
      int *active_individuals          = active_individuals_b.dev_pointer();

      auto q = (*this->scratch).get_queue();

      adadelta_krnl = std::make_unique<adadelta_kernel<queue_type>>(individuals_per_ligand,
                                                                  batch_ligands,
                                                                  batch_atoms,
                                                                  gradients_b,
                                                                  population_b,
                                                                  num_rotamers_p,
                                                                  adadelta_e_g2,
                                                                  adadelta_e_dw2,
                                                                  active_individuals,
                                                                  q,
                                                                  RHO,
                                                                  EPSILON);

      // Initialize the scoring kernel buffers
      this->score_stage->prepare(batch); //TODO L se non sbaglio l'ho aggiunto per quando deve fare solo local search nell'eseguibile stand alone
      geom_trans.prepare(batch);
    }

    void operator()() {
      assert(
          (((*this->scratch).template get<buffer_data_type::GRADIENTS>().num_elements() % batch_ligands) == 0) &&
          "Number of gradients is not a multiple of ligands in the batch");
      assert(adadelta_krnl && "Adadelta local search kernel method not yet prepared");

      auto &scores_b = (*this->scratch).template get<buffer_data_type::SCORES>();

      // TODO L try to move the reset of actives here which is more elegant, for now it is in adadelta cpp

      const bool only_local_search =
          ((*this->scratch).configuration.population_number == 1) &&
          ((*this->scratch).configuration.num_generations == 1);

      int j = 1;
      for (std::size_t i = 0; i < this->iterations; ++i) {
        geom_trans();
        
        if (only_local_search) {
          (*this->score_stage)();

          // copy scores back to host and print best score (first element)
          scores_b.copy_device2host();
          (*this->scratch).get_queue()->synchronize();
          if(i % (this->iterations/10) == 0){ // print every 10% of the process
            this->dump_pose(int(j++));
            printf("Iter: %ld, Score: %f\n", i, double(scores_b()[0]));
          }
        }

        // Fundamental part: compute gradient + adadelta update
        (this->score_stage).get()->compute_gradient();
        (*adadelta_krnl)();
      }

      if (only_local_search) {
        geom_trans();
        (*this->score_stage)();
      }
    }

    // TODO L Important: this was copied from genetic.hpp code, but not sure if it must be adapted 
    static std::size_t get_shared_ligand_mem(const int max_atoms, const knobs conf) {
      return 0;
    }

    static std::size_t get_private_ligand_mem(const int max_atoms, const knobs conf) {
      std::size_t mem{0};
      const int individuals_per_ligand = std::max(1, static_cast<int>(conf.population_number));

      mem += scoring_t<queue_type>::get_ligand_mem(max_atoms, conf);

      // One gradient per individual per ligand
      mem += sizeof(gradient) * individuals_per_ligand;
      // AdaDelta state buffers for each individual
      mem += sizeof(chromosome) * individuals_per_ligand; // E[g^2]
      mem += sizeof(chromosome) * individuals_per_ligand; // E[delta_w^2]
      mem += sizeof(int)        * individuals_per_ligand; // active flag
      return mem;
    }

    static int get_ligand_mem(const int max_atoms, const knobs conf) {
      return static_cast<int>(get_shared_ligand_mem(max_atoms, conf) +
                              get_private_ligand_mem(max_atoms, conf));
    }

    static batch_multiple get_batch_size(const int atoms,
                                         std::shared_ptr<queue_type> q,
                                         const knobs &conf,
                                         const size_t max_bucket_size) {
      (void) conf;
      const auto plain_multiple_info =
          normalize_batch_multiple(get_adadelta_batch_multiple<queue_type>(atoms, q));
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
    int batch_atoms;

    std::unique_ptr<adadelta_kernel<queue_type>> adadelta_krnl;
    geometric<queue_type> geom_trans;

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
