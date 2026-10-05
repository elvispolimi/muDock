/**
 * @file genetic_alpaka.cpp
 * @brief Implementation of the Genetic Algorithm (GA) kernels for Alpaka.
 * @details Implements the initialization (`initialize_alpaka`), reproduction/mutation (`iterate_alpaka`),
 *          and reduction (`finalize_alpaka`) kernels for evolving ligand conformations on accelerators.
 */

#include <alpaka/alpaka.hpp>
#include <cstdint>
#include <limits>
#include <mudock/alpaka_implementation/alpaka_random.hpp>
#include <mudock/alpaka_implementation/genetic_alpaka.hpp>
#include <mudock/alpaka_implementation/invoke_kernel_alpaka.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/utils.hpp>
#include <mudock/compute/devices_memory.hpp>
#include <mudock/alpaka_implementation/alpaka_random.hpp>

namespace mudock {
  /// @brief Thread-local cache for device-allocated random generator objects.
  thread_local device_memory<alpaka_random_object> alpaka_random_memory;

  namespace {
    static constexpr fp_type coordinate_step{static_cast<fp_type>(0.2)}; ///< Translation mutation scaling factor.
    static constexpr fp_type angle_step{4};                               ///< Rotation and torsion angle mutation scaling factor.

    /**
     * @brief Extracts a raw 32-bit random word from the generator state.
     * @param state Active PRNG engine state.
     * @return 32-bit unsigned integer.
     */
    ALPAKA_FN_ACC ALPAKA_FN_INLINE std::uint32_t next_random(alpaka_rand_state& state) {
      return static_cast<std::uint32_t>(state());
    }

    /**
     * @brief Generates a uniform floating-point value in [0, 1).
     * @details In debug mode (`is_debug()`), returns fixed constant `0.4f` for deterministic validation.
     * @param state Active PRNG engine state.
     * @return Uniform floating-point scalar.
     */
    ALPAKA_FN_ACC ALPAKA_FN_INLINE fp_type random_unit(alpaka_rand_state& state) {
      if constexpr (is_debug()) {
        return fp_type{0.4f};
      } else {
        return static_cast<fp_type>(next_random(state)) /
               static_cast<fp_type>(state.max());
      }
    }

    /**
     * @brief Generates a random value uniformly distributed between min and max.
     * @tparam T Numeric output type.
     * @param state Active PRNG engine state.
     * @param min Inclusive minimum bound.
     * @param max Inclusive maximum bound.
     * @return Uniformly distributed value of type T.
     */
    template<typename T>
    ALPAKA_FN_ACC ALPAKA_FN_INLINE T random_gen_alpaka(alpaka_rand_state& state, const T min, const T max) {
      return static_cast<T>((random_unit(state) * static_cast<fp_type>(max - min)) +
                            static_cast<fp_type>(min));
    }

    /// @brief Generates a random chromosome index for tournament participant selection.
    ALPAKA_FN_ACC ALPAKA_FN_INLINE int get_selection_distribution(alpaka_rand_state& state, const int population_number) {
      return random_gen_alpaka<int>(state, 0, population_number - 1);
    }

    /// @brief Generates initial conformation orientation perturbation in degrees [-45, 45].
    ALPAKA_FN_ACC ALPAKA_FN_INLINE fp_type get_init_change_distribution(alpaka_rand_state& state) {
      return random_gen_alpaka<fp_type>(state, -45, 45);
    }

    /// @brief Generates mutation perturbation step in [-10, 10].
    ALPAKA_FN_ACC ALPAKA_FN_INLINE fp_type get_mutation_change_distribution(alpaka_rand_state& state) {
      return random_gen_alpaka<fp_type>(state, -10, 10);
    }

    /// @brief Generates probability check value in [0, 1].
    ALPAKA_FN_ACC ALPAKA_FN_INLINE fp_type get_mutation_coin_distribution(alpaka_rand_state& state) {
      return random_gen_alpaka<fp_type>(state, 0, 1);
    }

    /// @brief Generates single-point crossover cut index across genes [0, 6 + num_rotamers].
    ALPAKA_FN_ACC ALPAKA_FN_INLINE int get_crossover_distribution(alpaka_rand_state& state, const int num_rotamers) {
      return random_gen_alpaka<int>(state, 0, 6 + num_rotamers);
    }

    /**
     * @brief Executes tournament selection to identify the fittest individual among candidates.
     * @param state Active PRNG state.
     * @param tournament_length Number of competing individuals in the tournament.
     * @param chromosome_number Total size of the candidate population.
     * @param[in] scores Energy score array for candidates.
     * @return Population index of the winning candidate with lowest energy.
     */
    ALPAKA_FN_ACC ALPAKA_FN_INLINE int tournament_selection_alpaka(alpaka_rand_state& state,
                                                  const int tournament_length,
                                                  const int chromosome_number,
                                                  const fp_type* __restrict__ scores) {
      int best_individual = get_selection_distribution(state, chromosome_number);
      for (int i = 0; i < tournament_length; ++i) {
        const auto contended = get_selection_distribution(state, chromosome_number);
        if (scores[contended] < scores[best_individual]) {
          best_individual = contended;
        }
      }
      return best_individual;
    }

    /**
     * @struct initialize_alpaka
     * @brief Device kernel generating initial random conformations for the population.
     * @details Populates chromosome gene values (translation and rotations) with initial
     *          random deviations and sets candidate scores to infinity.
     */
    struct initialize_alpaka {
      /**
       * @brief Kernel body initializing chromosomes and scores.
       * @tparam TAcc Alpaka accelerator type.
       * @param[in] acc Reference to the execution context.
       * @param chromosome_number Size of population per ligand.
       * @param[in] ligand_num_rotamers Rotatable bond counts per ligand.
       * @param[out] chromosomes Output array of generated chromosomes.
       * @param[out] ligand_scores Output array of candidate scores initialized to infinity.
       * @param[in,out] state Array of device PRNG states.
       */
      template<typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int chromosome_number,
                                    const int* __restrict__ ligand_num_rotamers,
                                    chromosome* __restrict__ chromosomes,
                                    fp_type* __restrict__ ligand_scores,
                                    alpaka_rand_state* __restrict__ state) const {
        const int ligand_id = static_cast<int>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        const int local_thread_id =
            static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        const int thread_per_block =
            static_cast<int>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);
        const int global_thread_id = local_thread_id + thread_per_block * ligand_id;

        const int num_rotamers = ligand_num_rotamers[ligand_id];
        chromosome* l_chromosomes = chromosomes + ligand_id * chromosome_number;
        fp_type* scores = ligand_scores + chromosome_number * ligand_id;
        
        alpaka_rand_local l_state = state[global_thread_id];

        for (int chromosome_index = local_thread_id; chromosome_index < chromosome_number;
             chromosome_index += thread_per_block) {
          scores[chromosome_index] = std::numeric_limits<fp_type>::infinity();

          chromosome& chromo = *(l_chromosomes + chromosome_index);
          ALPAKA_UNROLL(MUDOCK_UNROLL_FACTOR)
          for (int i{0}; i < 3; ++i) {
            chromo[i] = get_init_change_distribution(l_state) * coordinate_step;
          }
          ALPAKA_UNROLL(MUDOCK_UNROLL_FACTOR)
          for (int i{3}; i < 6 + num_rotamers; ++i) {
            chromo[i] = get_init_change_distribution(l_state) * angle_step;
          }
        }
        if constexpr (!std::is_reference_v<alpaka_rand_local>) {
          state[global_thread_id] = l_state;
        }
      }
    };

    /**
     * @struct iterate_alpaka
     * @brief Device kernel executing one generation of genetic search.
     * @details Performs tournament selection of parent pairs, single-point crossover,
     *          and stochastic coordinate/torsion mutation.
     */
    struct iterate_alpaka {
      /**
       * @brief Kernel body executing genetic selection, crossover, and mutation.
       * @tparam TAcc Alpaka accelerator type.
       * @param[in] acc Reference to the execution context.
       * @param tournament_length Size of tournament selection pool.
       * @param mutation_prob Probability of mutating each individual gene.
       * @param chromosome_number Size of population per ligand.
       * @param[in] ligand_num_rotamers Rotatable bond counts per ligand.
       * @param[in] chromosomes Current generation population.
       * @param[out] next_chromosomes Offspring population buffer.
       * @param[in] ligand_scores Fitness energy scores of current generation.
       * @param[in,out] state Array of device PRNG states.
       */
      template<typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int tournament_length,
                                    const fp_type mutation_prob,
                                    const int chromosome_number,
                                    const int* __restrict__ ligand_num_rotamers,
                                    chromosome* __restrict__ chromosomes,
                                    chromosome* __restrict__ next_chromosomes,
                                    fp_type* __restrict__ ligand_scores,
                                    alpaka_rand_state* __restrict__ state) const {
        const int ligand_id = static_cast<int>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        const int local_thread_id =
            static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        const int thread_per_block =
            static_cast<int>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);
        const int global_thread_id = local_thread_id + thread_per_block * ligand_id;

        const int num_rotamers = ligand_num_rotamers[ligand_id];
        chromosome* __restrict__ l_chromosomes = chromosomes + ligand_id * chromosome_number;
        chromosome* __restrict__ l_next_chromosomes = next_chromosomes + ligand_id * chromosome_number;
        const fp_type* __restrict__ scores = ligand_scores + chromosome_number * ligand_id;

        alpaka_rand_local l_state = state[global_thread_id];

        for (int chromosome_index = local_thread_id; chromosome_index < chromosome_number;
             chromosome_index += thread_per_block) {
          chromosome& next_chromosome = *(l_next_chromosomes + chromosome_index);

          const int best_individual_1 =
               tournament_selection_alpaka(l_state, tournament_length, chromosome_number, scores);
          const int best_individual_2 =
              tournament_selection_alpaka(l_state, tournament_length, chromosome_number, scores);

          const int split_index = get_crossover_distribution(l_state, num_rotamers);
          fp_type* dst = next_chromosome.data();
          const fp_type* p1 = l_chromosomes[best_individual_1].data();
          const fp_type* p2 = l_chromosomes[best_individual_2].data();
          for (int i = 0; i < (6 + num_rotamers); ++i) {
            dst[i] = (i < split_index) ? p1[i] : p2[i];
          }

          ALPAKA_UNROLL(MUDOCK_UNROLL_FACTOR)
          for (int i{0}; i < 3; ++i) {
            if (get_mutation_coin_distribution(l_state) < mutation_prob) {
              next_chromosome[i] += get_mutation_change_distribution(l_state) * coordinate_step;
            }
          }
          ALPAKA_UNROLL(MUDOCK_UNROLL_FACTOR)
          for (int i{3}; i < 6 + num_rotamers; ++i) {
            if (get_mutation_coin_distribution(l_state) < mutation_prob) {
              next_chromosome[i] += get_mutation_change_distribution(l_state) * angle_step;
            }
          }
        }
        if constexpr (!std::is_reference_v<alpaka_rand_local>) {
          state[global_thread_id] = l_state;
        }
      }
    };

    /**
     * @struct finalize_alpaka
     * @brief Device kernel locating the best-scoring candidate chromosome per ligand.
     * @details Employs parallel warp shuffle reduction (`alpaka::warp::shfl_down`) to
     *          find the minimum energy individual within each block and writes the result.
     */
    struct finalize_alpaka {
      /**
       * @brief Kernel body reducing candidate scores to find the minimum.
       * @tparam TAcc Alpaka accelerator type.
       * @param[in] acc Reference to the execution context.
       * @param chromosome_number Population size per ligand.
       * @param[in] ligand_num_rotamers Rotatable bond count per ligand.
       * @param[in] ligand_scores Array of evaluated scores.
       * @param[out] ligand_best_scores Array storing the lowest energy score per ligand.
       * @param[in] chromosomes Candidate chromosome population array.
       * @param[out] best_chromosomes Array storing the optimal winning chromosome per ligand.
       */
      template<typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int chromosome_number,
                                    const int* __restrict__ ligand_num_rotamers,
                                    fp_type* __restrict__ ligand_scores,
                                    fp_type* __restrict__ ligand_best_scores,
                                    chromosome* __restrict__ chromosomes,
                                    chromosome* __restrict__ best_chromosomes) const {
        const int ligand_id = static_cast<int>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        const int local_thread_id = static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        const int thread_per_block = static_cast<int>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);

        const int num_rotamers                 = ligand_num_rotamers[ligand_id];
        chromosome* __restrict__ l_chromosomes = chromosomes + ligand_id * chromosome_number;
        fp_type* __restrict__ scores           = ligand_scores + chromosome_number * ligand_id;

        int min_index = local_thread_id;
        fp_type min_score = min_index < chromosome_number ? scores[min_index] : std::numeric_limits<fp_type>::infinity();

        for (int chromosome_index = local_thread_id + thread_per_block; chromosome_index < chromosome_number;
             chromosome_index += thread_per_block) {
          if (min_score > scores[chromosome_index]) {
            min_index = chromosome_index;
            min_score = scores[chromosome_index];
          }
        }

        ALPAKA_UNROLL(MUDOCK_UNROLL_FACTOR)
        for (uint32_t offset = MUDOCK_ALPAKA_BLOCK_SIZE / 2; offset > 0; offset /= 2) {
          const fp_type other_min_score = alpaka::warp::shfl_down(acc, min_score, offset, MUDOCK_ALPAKA_BLOCK_SIZE);
          const int other_min_index     = alpaka::warp::shfl_down(acc, min_index, offset, MUDOCK_ALPAKA_BLOCK_SIZE);
          if (other_min_score < min_score) {
            min_score = other_min_score;
            min_index = other_min_index;
          }
        }

        if (local_thread_id == 0) {
          ligand_best_scores[ligand_id] = min_score;
          fp_type* dst                  = (best_chromosomes + ligand_id)->data();
          const fp_type* src            = (l_chromosomes + min_index)->data();
          for (int i = 0; i < (6 + num_rotamers); ++i) { dst[i] = src[i]; }
        }
      }
    };
  } // namespace

  /**
   * @brief Initializes population conformations on the Alpaka accelerator.
   * @details Allocates and seeds RNG state memory, then launches `initialize_alpaka`.
   */
  template<>
  void genetic_kernel<queue_alpaka>::initialize() {
    alpaka_random_memory.init(q);
    alpaka_random_memory.get_data()->alloc(batch_ligands * MUDOCK_ALPAKA_BLOCK_SIZE, seed);

    q->invoke_kernel<initialize_alpaka>(batch_ligands,
                                        population_number,
                                        num_rotamers_b,
                                        population,
                                        scores_b,
                                        alpaka_random_memory.get_data()->dev_pointer());
  }

  /**
   * @brief Advances candidate conformations by one genetic generation on the accelerator.
   * @details Launches `iterate_alpaka` to perform selection, crossover, and mutation in parallel.
   */
  template<>
  void genetic_kernel<queue_alpaka>::operator()() {
    q->invoke_kernel<iterate_alpaka>(batch_ligands,
                                     tournament_length,
                                     mutation_prob,
                                     population_number,
                                     num_rotamers_b,
                                     population,
                                     next_population,
                                     scores_b,
                                     alpaka_random_memory.get_data()->dev_pointer());
  }

  /**
   * @brief Selects the minimum energy candidate pose per ligand across the evaluated population.
   * @details Launches `finalize_alpaka` utilizing warp shuffle reduction.
   */
  template<>
  void genetic_kernel<queue_alpaka>::finalize() {
    q->invoke_kernel<finalize_alpaka>(batch_ligands,
                                      population_number,
                                      num_rotamers_b,
                                      scores_b,
                                      best_scores_b,
                                      population,
                                      best_chromosomes_b);
  }
} // namespace mudock
