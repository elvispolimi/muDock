#pragma once

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <fstream>
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

  // TODO L Important implement this
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

      this->standalone_local_search = ((*this->scratch).configuration.population_number == 1) &&
                                ((*this->scratch).configuration.num_generations == 1);
      
      // Initialize the scoring kernel buffers
      this->score_stage->prepare(batch); //TODO L se non sbaglio l'ho aggiunto per quando deve fare solo local search nell'eseguibile stand alone
      geom_trans.prepare(batch);
    }

    void operator()() {
      assert(
          (((*this->scratch).template get<buffer_data_type::GRADIENTS>().num_elements() % batch_ligands) == 0) &&
          "Number of gradients is not a multiple of ligands in the batch");
      assert(adadelta_krnl && "Adadelta local search kernel method not yet prepared");

      if (this->standalone_local_search) {
        run_standalone();
      } else {
        run_as_lga_step();
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
    int  batch_ligands;
    int  batch_atoms;
    std::unique_ptr<adadelta_kernel<queue_type>> adadelta_krnl;
    geometric<queue_type> geom_trans;

    fp_type get_ligand_com_displacement(std::optional<point3D>& initial_com) {
      assert(this->ligand_template.has_value() && "Ligand template not initialized for COM logging");

      auto& x_scratch_b = (*this->scratch).template get<buffer_data_type::X_SCRATCH>();
      auto& y_scratch_b = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>();
      auto& z_scratch_b = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>();
      x_scratch_b.copy_device2host();
      y_scratch_b.copy_device2host();
      z_scratch_b.copy_device2host();
      (*this->scratch).get_queue()->synchronize();

      const int num_atoms = this->ligand_template->num_atoms();
      const point3D current_com =
          compute_center_of_mass(x_scratch_b.host_pointer(), y_scratch_b.host_pointer(), z_scratch_b.host_pointer(), num_atoms);

      if (!initial_com) {
        initial_com = current_com;
        return fp_type{0};
      }

      return (current_com - *initial_com).magnitude();
    }

    void run_as_lga_step() {
      for (std::size_t i = 0; i < this->iterations; ++i) {
        geom_trans();
        (this->score_stage).get()->compute_gradient();
        (*adadelta_krnl)();
      }
    }

    void run_standalone() {
      auto &scores_b = (*this->scratch).template get<buffer_data_type::SCORES>();
      int dump_index = 1;

      const std::string score_log_path = "adadelta_scores.csv";
      const bool file_already_exists =
          std::filesystem::exists(score_log_path) && std::filesystem::file_size(score_log_path) > 0;
 
      std::ofstream score_log(score_log_path, std::ios::app);
      if (!file_already_exists) {
        score_log << "iteration,ligand,score\n";
      }

      const std::string com_log_path = "adadelta_com_distances.csv";
      const bool com_file_already_exists =
          std::filesystem::exists(com_log_path) && std::filesystem::file_size(com_log_path) > 0;
      std::ofstream com_log(com_log_path, std::ios::app);
      if (!com_file_already_exists) {
        com_log << "iteration,ligand,com_distance\n";
      }

      const std::string ligand_name = this->ligand_template ? this->ligand_template->properties.get(property_type::NAME) : std::string{"unknown"};
      
      const auto log_scores = [&](std::size_t iter) {
          score_log << iter << "," << ligand_name << "," << double(scores_b()[0]) << "\n";
      };

      std::optional<point3D> initial_com;
      const auto log_com_distance = [&](std::size_t iter) {
        const fp_type com_distance = get_ligand_com_displacement(initial_com);
        com_log << iter << "," << ligand_name << "," << double(com_distance) << "\n";
      };

      for (std::size_t iter = 0; iter < this->iterations; ++iter) {
        geom_trans();
        (*this->score_stage)();

        scores_b.copy_device2host();
        (*this->scratch).get_queue()->synchronize();

        log_scores(iter);

        this->dump_pose(dump_index++);

        log_com_distance(iter);

        (this->score_stage).get()->compute_gradient();
        (*adadelta_krnl)();
      }

      geom_trans();
      (*this->score_stage)();

      scores_b.copy_device2host();
      (*this->scratch).get_queue()->synchronize();
      log_scores(this->iterations);
      log_com_distance(this->iterations);
      score_log.close();
    }

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
