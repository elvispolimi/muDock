#pragma once

#include "mudock/type_alias.hpp"

#include <chrono>
#include <iostream>
#include <concepts>
#include <memory>
#include <mudock/chem/geom_ligand.hpp>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/log.hpp>
#if !defined(__CUDACC__) && !defined(__HIPCC__)
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scratchpad.hpp>
  #include <mudock/compute/convergence.hpp>
#endif
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cpp_implementation/mutate_cpp.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<typename queue_type>
  batch_multiple get_adadelta_conv_batch_multiple(const int, std::shared_ptr<queue_type>) {
    return {};
  }

  inline constexpr int MIN_ITERATIONS_TO_WARM_UP = 20;
  inline constexpr int MAX_BEST_SCORE_AGE = 5;

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct adadelta_convergence_kernel {
    static constexpr char adadelta_convergence_region_name[] = "adadelta_convergence";
    adadelta_convergence_kernel(const int batch_ligands_,
                                // const bool adadelta_autostop_,
                                const int population_number_,
                                int* __restrict__ converged_ligands_b_,
                                int* __restrict__ ls_selected_b_,
                                const fp_type *__restrict__ scores_b_,
                                fp_type* __restrict__ best_scores_b_,
                                int* __restrict__ best_scores_ages_b_,
                                std::shared_ptr<queue_type> q_)
                                : batch_ligands(batch_ligands_),
                                  // adadelta_autostop(adadelta_autostop_),
                                  population_number(population_number_),
                                  converged_ligands_b(converged_ligands_b_),
                                  ls_selected_b(ls_selected_b_),
                                  scores_b(scores_b_),
                                  best_scores_b(best_scores_b_),
                                  best_scores_ages_b(best_scores_ages_b_),
                                  q(q_) {}

    void operator()();

    adadelta_convergence_kernel(const adadelta_convergence_kernel&)            = default;
    adadelta_convergence_kernel(adadelta_convergence_kernel&&)                 = default;
    adadelta_convergence_kernel& operator=(const adadelta_convergence_kernel&) = delete;
    adadelta_convergence_kernel& operator=(adadelta_convergence_kernel&&)      = delete;

    ~adadelta_convergence_kernel() = default;

  private:
    const int batch_ligands;
    // const bool adadelta_autostop;
    const int population_number;
    int* __restrict__ converged_ligands_b;
    int* __restrict__ ls_selected_b;
    const fp_type *__restrict__ scores_b;
    fp_type* __restrict__ best_scores_b;
    int* __restrict__ best_scores_ages_b;
    std::shared_ptr<queue_type> q;
  };

#if !defined(__CUDACC__) && !defined(__HIPCC__)
  template<typename queue_t>
    requires std::derived_from<queue_t, queue>
  struct adadelta_convergence: public convergence<queue_t> {
    adadelta_convergence(std::shared_ptr<scratchpad<queue_t>> _scratch)
        : convergence<queue_t>(_scratch) {};

    void prepare(batch<static_molecule>& batch) {
      batch_ligands                = batch.num_ligands;
      auto q                       = (*this->scratch).get_queue();
      const auto population_number = std::max(1, static_cast<int>((*this->scratch).configuration.population_number));

      auto& converged_ligands_b = (*this->scratch).template get<buffer_data_type::CONVERGED_LIGANDS>();
      auto& ls_selected_b       = (*this->scratch).template get<buffer_data_type::LS_SELECTED>(); // recicle this buffer: must be reinitialized every generation
      auto& scores_b            = (*this->scratch).template get<buffer_data_type::SCORES>();
      auto& best_scores_b       = (*this->scratch).template get<buffer_data_type::BEST_SCORES>();
      auto& best_scores_ages_b  = (*this->scratch).template get<buffer_data_type::BEST_SCORES_AGE>();


      // initialize_buffer<buffer_data_type::CONVERGED_LIGANDS>(this->scratch, batch_ligands, 0, true);
      initialize_buffer<buffer_data_type::BEST_SCORES>(this->scratch, population_number * batch_ligands, big_bound<fp_type>(), true);
      initialize_buffer<buffer_data_type::BEST_SCORES_AGE>(this->scratch, population_number * batch_ligands, 0, true);

      int* __restrict__ converged_ligands_p = converged_ligands_b.dev_pointer();
      int* __restrict__ ls_selected_p       = ls_selected_b.dev_pointer();
      fp_type* __restrict__ scores_p        = scores_b.dev_pointer();
      fp_type* __restrict__ best_scores_p   = best_scores_b.dev_pointer();
      int* __restrict__ best_scores_ages_p  = best_scores_ages_b.dev_pointer();


      kernel = std::make_unique<adadelta_convergence_kernel<queue_t>>(batch_ligands,
                                                                    //  (*this->scratch).configuration.adadelta_autostop,
                                                                     population_number,
                                                                     converged_ligands_p,
                                                                     ls_selected_p,
                                                                     scores_p,
                                                                     best_scores_p,
                                                                     best_scores_ages_p,
                                                                     q);
    }

    void operator()() {
      // assert(
      //     ((*this->scratch).template get<buffer_data_type::CHROMOSOMES>().num_elements() % batch_ligands ==
      //      0) &&
      //     "Number of chromosomes per ligand is not a multiple of the number of ligands expected to be docked");

      assert(kernel && "Kernel method not yet prepared");
      (*kernel)();
    };


    // TODO fix function
    static std::size_t get_shared_ligand_mem(const int max_atoms, const knobs conf) {
      (void) max_atoms;
      (void) conf;
      return 0;
    }

    // TODO fix function
    static std::size_t get_private_ligand_mem(const int max_atoms, const knobs) {
      const int batch_rotamers              = max_atoms - 3;
      const int tot_rotamers_atoms_in_batch = max_atoms * batch_rotamers;
      std::size_t mem{0};
      return mem;
    }

    static int get_ligand_mem(const int max_atoms, const knobs conf) {
      return static_cast<int>(get_shared_ligand_mem(max_atoms, conf) +
                              get_private_ligand_mem(max_atoms, conf));
    }

    static batch_multiple get_batch_size(const int atoms,
                                         std::shared_ptr<queue_t> q,
                                         const knobs& conf,
                                         const size_t max_bucket_size) {
      (void) conf;
      const auto plain_multiple_info =
          normalize_batch_multiple(get_crystal_conv_batch_multiple<queue_t>(atoms, q));
      mudock::stage_bucket_trace("CONVERGENCE stage plain multiple for ",
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
    std::unique_ptr<adadelta_convergence_kernel<queue_t>> kernel;

    void teardown_impl(batch<static_molecule>& batch) {
      (void) batch;
    }
  };
#endif
} // namespace mudock
