#pragma once

#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <mudock/compute/queue.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/type_alias.hpp>
#if !defined(__CUDACC__) && !defined(__HIPCC__)
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/convergence.hpp>
#endif

namespace mudock {
  inline constexpr std::size_t MIN_ITERATIONS_TO_WARM_UP = 20;
  inline constexpr int MAX_BEST_SCORE_AGE = 5;

  // Floating-point predicates may be folded away under the Release build's
  // -ffast-math. Inspect the IEEE 754 exponent so invalid scores stay excluded.
#if defined(__CUDACC__) || defined(__HIPCC__)
  __host__ __device__
#endif
  inline bool finite_adadelta_score(const fp_type score) {
    static_assert(sizeof(fp_type) == sizeof(std::uint32_t) &&
                  std::numeric_limits<fp_type>::is_iec559);
    constexpr std::uint32_t exponent_mask = 0x7f800000U;
    return (std::bit_cast<std::uint32_t>(score) & exponent_mask) != exponent_mask;
  }

  // The score describes the current chromosome after completed_iterations updates.
  // Keep the best pose as well as its score: Adadelta need not improve monotonically.
#if defined(__CUDACC__) || defined(__HIPCC__)
  __host__ __device__
#endif
  inline void update_adadelta_convergence(const fp_type score,
                                          fp_type& best_score,
                                          int& best_score_age,
                                          int& ls_selected,
                                          chromosome& current,
                                          chromosome& best,
                                          const std::size_t completed_iterations,
                                          const bool finalize) {
    if (!ls_selected) return;

    if (finite_adadelta_score(score) && score < best_score) {
      best_score = score;
      best_score_age = 0;
      best = current;
    } else if (best_score_age < MAX_BEST_SCORE_AGE) {
      ++best_score_age;
    }

    if (finalize || (completed_iterations >= MIN_ITERATIONS_TO_WARM_UP &&
                     best_score_age >= MAX_BEST_SCORE_AGE)) {
      if (finite_adadelta_score(best_score)) current = best;
      ls_selected = 0;
    }
  }

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct adadelta_convergence_kernel {
    static constexpr char adadelta_convergence_region_name[] = "adadelta_convergence";

    adadelta_convergence_kernel(const int batch_ligands_,
                                const int population_number_,
                                const int* converged_ligands_b_,
                                int* ls_selected_b_,
                                const fp_type* scores_b_,
                                fp_type* best_scores_b_,
                                int* best_scores_ages_b_,
                                chromosome* chromosomes_b_,
                                chromosome* best_chromosomes_b_,
                                std::shared_ptr<queue_type> q_)
        : batch_ligands(batch_ligands_),
          population_number(population_number_),
          converged_ligands_b(converged_ligands_b_),
          ls_selected_b(ls_selected_b_),
          scores_b(scores_b_),
          best_scores_b(best_scores_b_),
          best_scores_ages_b(best_scores_ages_b_),
          chromosomes_b(chromosomes_b_),
          best_chromosomes_b(best_chromosomes_b_),
          q(q_) {}

    void operator()(std::size_t completed_iterations, bool finalize);
    void set_chromosomes_buffer(chromosome* chromosomes) { chromosomes_b = chromosomes; }

  private:
    const int batch_ligands;
    const int population_number;
    const int* converged_ligands_b;
    int* ls_selected_b;
    const fp_type* scores_b;
    fp_type* best_scores_b;
    int* best_scores_ages_b;
    chromosome* chromosomes_b;
    chromosome* best_chromosomes_b;
    std::shared_ptr<queue_type> q;
  };

#if !defined(__CUDACC__) && !defined(__HIPCC__)
  template<typename queue_t>
    requires std::derived_from<queue_t, queue>
  struct adadelta_convergence: public convergence<queue_t> {
    explicit adadelta_convergence(std::shared_ptr<scratchpad<queue_t>> scratch)
        : convergence<queue_t>(scratch) {}

    void prepare(batch<static_molecule>& batch) override {
      const auto population_number = std::max(1, static_cast<int>(this->scratch->configuration.population_number));
      individual_count = static_cast<std::size_t>(population_number) * static_cast<std::size_t>(batch.num_ligands);

      auto& converged = this->scratch->template get<buffer_data_type::CONVERGED_LIGANDS>();
      auto& selected = this->scratch->template get<buffer_data_type::LS_SELECTED>();
      auto& scores = this->scratch->template get<buffer_data_type::SCORES>();
      auto& best_scores = this->scratch->template get<buffer_data_type::BEST_SCORES>();
      auto& ages = this->scratch->template get<buffer_data_type::BEST_SCORES_AGE>();
      auto& chromosomes = this->scratch->template get<buffer_data_type::CHROMOSOMES>();
      auto& best_chromosomes = this->scratch->template get<buffer_data_type::ADADELTA_BEST_CHROMOSOMES>();

      best_chromosomes.alloc(individual_count);
      reset();
      kernel = std::make_unique<adadelta_convergence_kernel<queue_t>>(batch.num_ligands,
                                                                    population_number,
                                                                    converged.dev_pointer(),
                                                                    selected.dev_pointer(),
                                                                    scores.dev_pointer(),
                                                                    best_scores.dev_pointer(),
                                                                    ages.dev_pointer(),
                                                                    chromosomes.dev_pointer(),
                                                                    best_chromosomes.dev_pointer(),
                                                                    this->scratch->get_queue());
    }

    // Population slots hold new individuals each generation. Their convergence
    // history must start afresh along with the Adadelta accumulators and mask.
    void reset() {
      completed_iterations = 0;
      initialize_buffer<buffer_data_type::BEST_SCORES>(this->scratch, individual_count,
                                                       std::numeric_limits<fp_type>::infinity(), true);
      initialize_buffer<buffer_data_type::BEST_SCORES_AGE>(this->scratch, individual_count, 0, true);
    }

    void operator()() override {
      assert(kernel && "Adadelta convergence kernel not yet prepared");
      (*kernel)(completed_iterations++, false);
    }

    // Evaluate the last update before restoring the best pose at the budget limit.
    void finalize() {
      assert(kernel && "Adadelta convergence kernel not yet prepared");
      (*kernel)(completed_iterations, true);
    }

    void set_chromosomes_buffer(chromosome* chromosomes) {
      if (kernel) kernel->set_chromosomes_buffer(chromosomes);
    }

    static std::size_t get_private_ligand_mem(const int, const knobs conf) {
      return static_cast<std::size_t>(std::max(1, static_cast<int>(conf.population_number))) *
             (sizeof(fp_type) + sizeof(int) + sizeof(chromosome));
    }

  private:
    std::size_t individual_count{0};
    std::size_t completed_iterations{0};
    std::unique_ptr<adadelta_convergence_kernel<queue_t>> kernel;

    void teardown_impl(batch<static_molecule>&) override {}
  };
#endif
} // namespace mudock
