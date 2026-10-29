#pragma once

#include <vector>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <algorithm>
#include <numeric>
#include <mudock/batch.hpp>
#include <mudock/csv_logger.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/compute/adadelta_kernel.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/compute/ligand_placement.hpp>
#if !defined(__CUDACC__) && !defined(__HIPCC__)
  #include <mudock/format/reader.hpp>
  #include <mudock/compute/local_search.hpp>
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>

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
    // static constexpr fp_type RHO     = 0.85f;
    // static constexpr fp_type EPSILON = 1e-6f;
    
    adadelta(std::shared_ptr<scratchpad<queue_type>> _scratch,
             std::shared_ptr<scoring_t<queue_type>> _score,
             ligand_placement _placement = {})
             : local_search<queue_type, scoring_t>(_scratch, _score),
               geom_trans(_scratch, _score->get_protein(), _placement) {}

    void prepare(batch<static_molecule> &batch) {
      const knobs& configuration = (*this->scratch).configuration;
      // TODO L check if this makes a copy or a reference
      if (configuration.crystal_path) {
        this->crystal_template = mudock::parser<static_molecule>(std::filesystem::path(*configuration.crystal_path));
      } else {
        this->crystal_template = *batch.molecules[0];
      }

      if (this->crystal_template->num_atoms() != batch.molecules[0]->num_atoms()) {
        throw std::runtime_error("Crystal and docked ligand have different atom counts");
      }
      batch_ligands = batch.num_ligands;
      batch_atoms   = batch.batch_max_atoms;
      individuals_per_ligand = std::max(1, static_cast<int>(configuration.population_number));
      local_search_rate = static_cast<int>(configuration.lsrate);
      local_search_on_best = configuration.ls_on_best;
      rho = configuration.rho;
      epsilon = configuration.epsilon;
      this->standalone_local_search = this->is_standalone_local_search(configuration);
      initial_scores.resize(batch_ligands);
      initial_rmsds.resize(batch_ligands);
      final_scores.resize(batch_ligands);
      final_rmsds.resize(batch_ligands);
      this->ligand_names.resize(batch_ligands);
      for (int ligand_index = 0; ligand_index < batch_ligands; ++ligand_index) {
        this->ligand_names[ligand_index] = batch.molecules[ligand_index]->properties.get(property_type::NAME);
      }


      auto &adadelta_e_g2_b  = (*this->scratch).template get<buffer_data_type::ADADELTA_E_G2>();
      auto &adadelta_e_dw2_b = (*this->scratch).template get<buffer_data_type::ADADELTA_E_DW2>();
      auto &chromosomes_b    = (*this->scratch).template get<buffer_data_type::CHROMOSOMES>();
      auto &gradients_b       = (*this->scratch).template get<buffer_data_type::GRADIENTS>();
      const size_t gradient_count = static_cast<size_t>(batch_ligands) * static_cast<size_t>(individuals_per_ligand);


      auto& num_rotamers_b = (*this->scratch).template get<buffer_data_type::NUM_ROTAMERS>();
      num_rotamers_b.alloc(batch_ligands);
      load_num_rotamers<queue_type>(batch, this->scratch);
      
      // AdaDelta state: needed in both modes, and must always start from zero
      initialize_buffer<buffer_data_type::ADADELTA_E_G2>(this->scratch, gradient_count, chromosome{}, true);
      initialize_buffer<buffer_data_type::ADADELTA_E_DW2>(this->scratch, gradient_count, chromosome{}, true);
      initialize_buffer<buffer_data_type::GRADIENTS>(this->scratch, gradient_count, gradient{}, true);

      if (this->standalone_local_search) {
        initialize_buffer<buffer_data_type::CHROMOSOMES>(this->scratch, gradient_count, chromosome{}, true);
      }
      
      // TODO L move active from lsrate initialization from here to local search generic or LGA?
      // Initialize active population
      auto &active_individuals_b = (*this->scratch).template get<buffer_data_type::ACTIVE_INDIVIDUALS>();
      if (!active_individuals_b.is_valid() || active_individuals_b.num_elements() != gradient_count) {
        active_individuals_b.alloc(gradient_count);
        active_individuals_b.set_valid();
      }
      std::vector<int> active_init(gradient_count);
      
      std::mt19937 rng(configuration.seed.value_or(std::random_device{}()));
      std::bernoulli_distribution dist(static_cast<double>(static_cast<fp_type>(local_search_rate) / fp_type{100}));

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
      chromosome *chromosomes          = chromosomes_b.dev_pointer();
      gradient   *gradients            = gradients_b.dev_pointer();
      int        *active_individuals   = active_individuals_b.dev_pointer();

      auto q = (*this->scratch).get_queue();

      auto &converged_ligands_b = (*this->scratch).template get<buffer_data_type::CONVERGED_LIGANDS>();
      if (!converged_ligands_b.is_valid() || converged_ligands_b.num_elements() != static_cast<size_t>(batch_ligands)) {
        converged_ligands_b.alloc(batch_ligands);
        converged_ligands_b.set_valid();
      }
      converged_ligands_b.copy_host2device();
      int *converged_ligands = converged_ligands_b.dev_pointer();

      adadelta_krnl = std::make_unique<adadelta_kernel<queue_type>>(individuals_per_ligand,
                                                                    batch_ligands,
                                                                    batch_atoms,
                                                                    gradients,
                                                                    chromosomes,
                                                                    num_rotamers_p,
                                                                    adadelta_e_g2,
                                                                    adadelta_e_dw2,
                                                                    active_individuals,
                                                                    converged_ligands,
                                                                    q,
                                                                    rho,
                                                                    epsilon);

      
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

    inline void set_chromosomes_buffer(chromosome* chromosomes_b_) {
      if (adadelta_krnl) {
        adadelta_krnl->set_chromosomes_buffer(chromosomes_b_);
        geom_trans.set_chromosomes_buffer(chromosomes_b_);
      }
    }

    static std::size_t get_shared_ligand_mem(const int max_atoms, const knobs conf) {
      (void) max_atoms;
      (void) conf;
      return 0;
    }

    static std::size_t get_private_ligand_mem(const int max_atoms, const knobs conf) {
      std::size_t mem{0};

      const int individuals_per_ligand =
          std::max(1, static_cast<int>(conf.population_number));

      // Memory required by the scoring stage for one ligand.
      mem += scoring_t<queue_type>::get_ligand_mem(max_atoms, conf);

      // One gradient per individual.
      mem += sizeof(gradient) * individuals_per_ligand;

      // One AdaDelta E[g^2] state per individual.
      mem += sizeof(chromosome) * individuals_per_ligand;

      // One AdaDelta E[delta_w^2] state per individual.
      mem += sizeof(chromosome) * individuals_per_ligand;

      // One active flag per individual.
      mem += sizeof(int) * individuals_per_ligand;

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
    int individuals_per_ligand;
    int local_search_rate;
    bool local_search_on_best;
    fp_type rho;
    fp_type epsilon;
    std::vector<fp_type> initial_scores;
    std::vector<fp_type> initial_rmsds;
    std::vector<fp_type> final_scores;
    std::vector<fp_type> final_rmsds;

    std::unique_ptr<adadelta_kernel<queue_type>> adadelta_krnl;
    geometric<queue_type> geom_trans;

    // TODO L valid only for standalone (one individual)
    fp_type get_ligand_rmsd(int ligand_index, const fp_type* x_scratch_b, const fp_type* y_scratch_b, const fp_type* z_scratch_b) {
      assert(this->crystal_template.has_value() && "Ligand template not initialized for RMSD logging");

      const int num_atoms   = this->crystal_template->num_atoms();
      const int atom_stride = ligand_index * batch_atoms;

      const fp_type *__restrict__ scratch_x = x_scratch_b + atom_stride;
      const fp_type *__restrict__ scratch_y = y_scratch_b + atom_stride;
      const fp_type *__restrict__ scratch_z = z_scratch_b + atom_stride;


      return compute_rmsd((*this->crystal_template).x(), (*this->crystal_template).y(), (*this->crystal_template).z(),
                          scratch_x, scratch_y, scratch_z, 
                          num_atoms);
    }



    void run_as_lga_step() {

      // Reinitialize buffers
      const size_t gradient_count = static_cast<size_t>(batch_ligands) * static_cast<size_t>(individuals_per_ligand);
      initialize_buffer<buffer_data_type::ADADELTA_E_G2>(this->scratch, gradient_count, chromosome{}, true);
      initialize_buffer<buffer_data_type::ADADELTA_E_DW2>(this->scratch, gradient_count, chromosome{}, true);

      if (local_search_on_best) {
        // Mark top n individuals based on lsrate
        auto& scores_b             = (*this->scratch).template get<buffer_data_type::SCORES>();
        auto& active_individuals_b = (*this->scratch).template get<buffer_data_type::ACTIVE_INDIVIDUALS>();
        scores_b.copy_device2host();
        active_individuals_b.copy_device2host();
        (*this->scratch).get_queue()->synchronize();
        fp_type* __restrict__ scores_p  = scores_b.host_pointer();
        int* __restrict__ active_individuals_p  = active_individuals_b.host_pointer();
        for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
          int     *__restrict__ active_individuals_l = active_individuals_p + ligand_index * individuals_per_ligand;
          fp_type *__restrict__ scores_l             = scores_p + ligand_index * individuals_per_ligand;
          const int n = static_cast<int>(static_cast<fp_type>(local_search_rate) * static_cast<fp_type>(individuals_per_ligand) / fp_type{100});
          markTopNActive(scores_l, active_individuals_l, n, individuals_per_ligand);
        }
      }

      // Run local search
      for (std::size_t i = 0; i < this->iterations; ++i) {
        geom_trans();
        (this->score_stage).get()->compute_gradient();
        (*adadelta_krnl)();
      }

      ///////////////////////////////////////////////////////////////////////////////////////////
      // TODO L : Experiments on rmsd, not needed in general for computation of local search
      // WARNING: GETTING SCORES ASSUMING THERE IS ONLY ONE LIGAND IN THE BATCH
      fp_type* __restrict__ scores = (*this->scratch).template get<buffer_data_type::SCORES>().dev_pointer();
      int best_index = 0;
      for(int i = 0; i < individuals_per_ligand; ++i){
        if (scores[i] < scores[best_index]){
          best_index = i;
        }
      }
      // const std::string ligand_name = this->crystal_template ? this->crystal_template->properties.get(property_type::NAME) : std::string{"unknown"};
      // remove this code after experiments ^
      ///////////////////////////////////////////////////////////////////////////////////////////

    }

    // scores: input array of length m
    // active_individuals: output array of length m (must be pre-allocated by caller)
    // m: length of both arrays
    // n: how many top scores to mark active
    template <typename fp_type>
    void markTopNActive(const fp_type* scores, int* active_individuals, int n, int m) {
        std::fill(active_individuals, active_individuals + m, 0);
        if (n <= 0) return;
        if (n >= m) { std::fill(active_individuals, active_individuals + m, 1); return; }

        // index array on the heap (or use alloca/stack array if m is small/fixed)
        int* idx = new int[m];
        std::iota(idx, idx + m, 0);  // idx = [0, 1, 2, ..., m-1]

        std::nth_element(idx, idx + n, idx + m,
            [scores](int a, int b) { return scores[a] < scores[b]; });

        for (int i = 0; i < n; i++) {
            active_individuals[idx[i]] = 1;
        }

        delete[] idx;
    }

    void run_standalone() {

      // What to log during the local search
      // csv_logger score_logger("adadelta_scores.csv", {"ligand", "iteration", "standalone_ls_score"});
      // csv_logger rmsd_logger("adadelta_rmsd.csv", {"ligand", "iteration", "rmsd"});
      // std::string complex_id = this->crystal_template ? this->crystal_template->properties.get(property_type::NAME) : std::string{"unknown"};
      // const auto pos = complex_id.find('_');
      // if (pos != std::string::npos) complex_id.erase(pos);
      // static csv_logger crystal_experiment_logger(complex_id + ".csv", {"ligand_index", "initial_score", "final_score", "initial_rmsd", "final_rmsd", "iterations"});
      
      // auto &scores_b = (*this->scratch).template get<buffer_data_type::SCORES>();
      // fp_type crystal_score = 0.0;

      // WARNING: THIS IS VALID ONLY FOR ONE LIGAND AT A TIME
      // const int num_atoms = this->crystal_template->num_atoms();
      // const long int num_rotamers = this->crystal_template->num_rotamers();
      

      for (std::size_t iter = 1; iter <= this->iterations; ++iter) {
        geom_trans();
        (*this->score_stage)();
        (this->score_stage).get()->compute_gradient();
        (*adadelta_krnl)();

        // Dump pose for visualization.
        // this->dump_pose(iter);
        // if (iter == 1) {
        //   auto& scores_b = (*this->scratch).template get<buffer_data_type::SCORES>();
        //   auto& x_scratch_b = (*this->scratch).template get<buffer_data_type::X_SCRATCH>();
        //   auto& y_scratch_b = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>();
        //   auto& z_scratch_b = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>();
        //   scores_b.copy_device2host();
        //   x_scratch_b.copy_device2host();
        //   y_scratch_b.copy_device2host();
        //   z_scratch_b.copy_device2host();
        //   (*this->scratch).get_queue()->synchronize();
        //   for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        //     initial_scores[ligand_index] = scores_b()[ligand_index];
        //     initial_rmsds[ligand_index] = get_ligand_rmsd(ligand_index, x_scratch_b.host_pointer(), y_scratch_b.host_pointer(), z_scratch_b.host_pointer());
        //   }
        // }
          
      }
      
      geom_trans();
      (*this->score_stage)();
      
      // auto& scores_b = (*this->scratch).template get<buffer_data_type::SCORES>();
      // auto& x_scratch_b = (*this->scratch).template get<buffer_data_type::X_SCRATCH>();
      // auto& y_scratch_b = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>();
      // auto& z_scratch_b = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>();
      // scores_b.copy_device2host();
      // x_scratch_b.copy_device2host();
      // y_scratch_b.copy_device2host();
      // z_scratch_b.copy_device2host();
      // const fp_type* x_scratch_b_ = x_scratch_b.host_pointer();
      // const fp_type* y_scratch_b_ = y_scratch_b.host_pointer();
      // const fp_type* z_scratch_b_ = z_scratch_b.host_pointer();

      // (*this->scratch).get_queue()->synchronize();
      // for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      //   final_scores[ligand_index] = scores_b()[ligand_index];
      //   final_rmsds[ligand_index] = get_ligand_rmsd(ligand_index, x_scratch_b_, y_scratch_b_, z_scratch_b_);
      //   crystal_experiment_logger.log(this->ligand_names[ligand_index], initial_scores[ligand_index], final_scores[ligand_index], initial_rmsds[ligand_index], final_rmsds[ligand_index], this->iterations);
      //   this->dump_pose(complex_id, ligand_index, x_scratch_b_, y_scratch_b_, z_scratch_b_, batch_atoms);
      // }

      

    }

    // TODO L: Implement teardown
    void teardown_impl(batch<static_molecule> &batch) override {
      if (this->standalone_local_search){
        auto &scores_b              = (*this->scratch).template get<buffer_data_type::SCORES>();
        const int scores_per_ligand = static_cast<int>(scores_b.num_elements() / batch_ligands);
        scores_b.copy_device2host();
        (*this->scratch).get_queue()->synchronize();
        for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
          auto &ligand          = *batch.molecules[ligand_index];
          // const int num_atoms = ligand.num_atoms();
          // const long int num_rotamers = ligand.num_rotamers();
          const int score_index = ligand_index * scores_per_ligand;
          ligand.properties.assign(property_type::SCORE, std::to_string(scores_b()[score_index]));
          // ligand.properties.assign(property_type::NUM_ATOMS, std::to_string(num_atoms));
          // ligand.properties.assign(property_type::NUM_ROT, std::to_string(num_rotamers));
          ligand.properties.assign(property_type::ITER, std::to_string(this->iterations));
          ligand.properties.assign(property_type::RMSD_BEST_SCORING_POSE, std::to_string(final_rmsds[ligand_index]));
        }
      }
    };

  };
  #endif
} // namespace mudock
