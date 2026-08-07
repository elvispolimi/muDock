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
  thread_local device_memory<alpaka_random_object> alpaka_random_memory;

  namespace {
    static constexpr fp_type coordinate_step{static_cast<fp_type>(0.2)};
    static constexpr fp_type angle_step{4};

    ALPAKA_FN_ACC ALPAKA_FN_INLINE std::uint32_t next_random(alpaka_rand_state& state) {
      return state();
    }

    ALPAKA_FN_ACC ALPAKA_FN_INLINE fp_type random_unit(alpaka_rand_state& state) {
      if constexpr (is_debug()) {
        return fp_type{0.4f};
      } else {
        return static_cast<fp_type>(next_random(state)) /
               static_cast<fp_type>(state.max());
      }
    }

    template<typename T>
    ALPAKA_FN_ACC ALPAKA_FN_INLINE T random_gen_alpaka(alpaka_rand_state& state, const T min, const T max) {
      return static_cast<T>((random_unit(state) * static_cast<fp_type>(max - min)) +
                            static_cast<fp_type>(min));
    }

    ALPAKA_FN_ACC ALPAKA_FN_INLINE int get_selection_distribution(alpaka_rand_state& state, const int population_number) {
      return random_gen_alpaka<int>(state, 0, population_number - 1);
    }

    ALPAKA_FN_ACC ALPAKA_FN_INLINE fp_type get_init_change_distribution(alpaka_rand_state& state) {
      return random_gen_alpaka<fp_type>(state, -45, 45);
    }

    ALPAKA_FN_ACC ALPAKA_FN_INLINE fp_type get_mutation_change_distribution(alpaka_rand_state& state) {
      return random_gen_alpaka<fp_type>(state, -10, 10);
    }

    ALPAKA_FN_ACC ALPAKA_FN_INLINE fp_type get_mutation_coin_distribution(alpaka_rand_state& state) {
      return random_gen_alpaka<fp_type>(state, 0, 1);
    }

    ALPAKA_FN_ACC ALPAKA_FN_INLINE int get_crossover_distribution(alpaka_rand_state& state, const int num_rotamers) {
      return random_gen_alpaka<int>(state, 0, 6 + num_rotamers);
    }

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

    struct initialize_alpaka {
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
        
        alpaka_rand_state l_state = state[global_thread_id];

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
        state[global_thread_id] = l_state;
      }
    };

    struct iterate_alpaka {
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

        alpaka_rand_state l_state = state[global_thread_id];

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
        state[global_thread_id] = l_state;
      }
    };

    struct finalize_alpaka {
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

  template<>
  void genetic_kernel<queue_alpaka>::initialize() {
    alpaka_random_memory.init(q);
    alpaka_random_memory.get_data()->alloc(batch_ligands * MUDOCK_ALPAKA_BLOCK_SIZE, seed);

    q->invoke_kernel<initialize_alpaka>(batch_ligands,
                                        MUDOCK_ALPAKA_BLOCK_SIZE,
                                        population_number,
                                        num_rotamers_b,
                                        population,
                                        scores_b,
                                        alpaka_random_memory.get_data()->dev_pointer());
  }

  template<>
  void genetic_kernel<queue_alpaka>::operator()() {
    q->invoke_kernel<iterate_alpaka>(batch_ligands,
                                     MUDOCK_ALPAKA_BLOCK_SIZE,
                                     tournament_length,
                                     mutation_prob,
                                     population_number,
                                     num_rotamers_b,
                                     population,
                                     next_population,
                                     scores_b,
                                     alpaka_random_memory.get_data()->dev_pointer());
  }

  template<>
  void genetic_kernel<queue_alpaka>::finalize() {
    q->invoke_kernel<finalize_alpaka>(batch_ligands,
                                      MUDOCK_ALPAKA_BLOCK_SIZE,
                                      population_number,
                                      num_rotamers_b,
                                      scores_b,
                                      best_scores_b,
                                      population,
                                      best_chromosomes_b);
  }
} // namespace mudock
