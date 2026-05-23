#include <mudock/alpaka_implementation/genetic_alpaka.hpp>
#include <mudock/alpaka_implementation/invoke_kernel_alpaka.hpp>

#include <alpaka/alpaka.hpp>

#include <cstdint>
#include <limits>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/utils.hpp>

#ifndef MUDOCK_ALPAKA_BLOCK_SIZE
  #define MUDOCK_ALPAKA_BLOCK_SIZE 32
#endif

namespace mudock {
  namespace {
    static constexpr fp_type coordinate_step{static_cast<fp_type>(0.2)};
    static constexpr fp_type angle_step{4};

    ALPAKA_FN_ACC std::uint32_t next_random(std::uint32_t& state) {
      state ^= state << 13;
      state ^= state >> 17;
      state ^= state << 5;
      return state;
    }

    ALPAKA_FN_ACC fp_type random_unit(std::uint32_t& state) {
      if constexpr (is_debug()) {
        return fp_type{0.4};
      } else {
        return static_cast<fp_type>(next_random(state)) /
               static_cast<fp_type>(std::numeric_limits<std::uint32_t>::max());
      }
    }

    template<typename T>
    ALPAKA_FN_ACC T random_gen_alpaka(std::uint32_t& state, const T min, const T max) {
      return static_cast<T>((random_unit(state) * static_cast<fp_type>(max - min)) +
                            static_cast<fp_type>(min));
    }

    ALPAKA_FN_ACC int get_selection_distribution(std::uint32_t& state, const int population_number) {
      return random_gen_alpaka<int>(state, 0, population_number - 1);
    }

    ALPAKA_FN_ACC fp_type get_init_change_distribution(std::uint32_t& state) {
      return random_gen_alpaka<fp_type>(state, -45, 45);
    }

    ALPAKA_FN_ACC fp_type get_mutation_change_distribution(std::uint32_t& state) {
      return random_gen_alpaka<fp_type>(state, -10, 10);
    }

    ALPAKA_FN_ACC fp_type get_mutation_coin_distribution(std::uint32_t& state) {
      return random_gen_alpaka<fp_type>(state, 0, 1);
    }

    ALPAKA_FN_ACC int get_crossover_distribution(std::uint32_t& state, const int num_rotamers) {
      return random_gen_alpaka<int>(state, 0, 6 + num_rotamers);
    }

    ALPAKA_FN_ACC int tournament_selection_alpaka(std::uint32_t& state,
                                                  const int tournament_length,
                                                  const int chromosome_number,
                                                  const fp_type* scores) {
      int best_individual = get_selection_distribution(state, chromosome_number);
      for (int i = 0; i < tournament_length; ++i) {
        const auto contended = get_selection_distribution(state, chromosome_number);
        if (scores[contended] < scores[best_individual]) {
          best_individual = contended;
        }
      }
      return best_individual;
    }

    ALPAKA_FN_ACC std::uint32_t make_seed(const std::size_t seed,
                                          const int ligand_id,
                                          const int local_thread_id,
                                          const int salt) {
      auto state = static_cast<std::uint32_t>(seed) ^ 0x9e3779b9u;
      state ^= static_cast<std::uint32_t>(ligand_id + 1) * 0x85ebca6bu;
      state ^= static_cast<std::uint32_t>(local_thread_id + 1) * 0xc2b2ae35u;
      state ^= static_cast<std::uint32_t>(salt + 1) * 0x27d4eb2fu;
      return state ? state : 0x6d2b79f5u;
    }

    struct initialize_alpaka {
      template<typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int chromosome_number,
                                    const int* ligand_num_rotamers,
                                    chromosome* chromosomes,
                                    fp_type* ligand_scores,
                                    const std::size_t seed) const {
        const int ligand_id = static_cast<int>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        const int local_thread_id =
            static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        const int thread_per_block =
            static_cast<int>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);

        const int num_rotamers = ligand_num_rotamers[ligand_id];
        chromosome* l_chromosomes = chromosomes + ligand_id * chromosome_number;
        fp_type* scores = ligand_scores + chromosome_number * ligand_id;

        for (int chromosome_index = local_thread_id; chromosome_index < chromosome_number;
             chromosome_index += thread_per_block) {
          scores[chromosome_index] = std::numeric_limits<fp_type>::infinity();

          auto state = make_seed(seed, ligand_id, local_thread_id, chromosome_index);
          chromosome& chromo = *(l_chromosomes + chromosome_index);
          for (int i{0}; i < 3; ++i) {
            chromo[i] = get_init_change_distribution(state) * coordinate_step;
          }
          for (int i{3}; i < 6 + num_rotamers; ++i) {
            chromo[i] = get_init_change_distribution(state) * angle_step;
          }
        }
      }
    };

    struct iterate_alpaka {
      template<typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int tournament_length,
                                    const fp_type mutation_prob,
                                    const int chromosome_number,
                                    const int* ligand_num_rotamers,
                                    chromosome* chromosomes,
                                    chromosome* next_chromosomes,
                                    fp_type* ligand_scores,
                                    const std::size_t seed) const {
        const int ligand_id = static_cast<int>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        const int local_thread_id =
            static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        const int thread_per_block =
            static_cast<int>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);

        const int num_rotamers = ligand_num_rotamers[ligand_id];
        chromosome* l_chromosomes = chromosomes + ligand_id * chromosome_number;
        chromosome* l_next_chromosomes = next_chromosomes + ligand_id * chromosome_number;
        const fp_type* scores = ligand_scores + chromosome_number * ligand_id;

        for (int chromosome_index = local_thread_id; chromosome_index < chromosome_number;
             chromosome_index += thread_per_block) {
          auto state = make_seed(seed, ligand_id, local_thread_id, chromosome_index);
          chromosome& next_chromosome = *(l_next_chromosomes + chromosome_index);

          const int best_individual_1 =
              tournament_selection_alpaka(state, tournament_length, chromosome_number, scores);
          const int best_individual_2 =
              tournament_selection_alpaka(state, tournament_length, chromosome_number, scores);

          const int split_index = get_crossover_distribution(state, num_rotamers);
          fp_type* dst = next_chromosome.data();
          const fp_type* p1 = l_chromosomes[best_individual_1].data();
          const fp_type* p2 = l_chromosomes[best_individual_2].data();
          for (int i = 0; i < (6 + num_rotamers); ++i) {
            dst[i] = (i < split_index) ? p1[i] : p2[i];
          }

          for (int i{0}; i < 3; ++i) {
            if (get_mutation_coin_distribution(state) < mutation_prob) {
              next_chromosome[i] += get_mutation_change_distribution(state) * coordinate_step;
            }
          }
          for (int i{3}; i < 6 + num_rotamers; ++i) {
            if (get_mutation_coin_distribution(state) < mutation_prob) {
              next_chromosome[i] += get_mutation_change_distribution(state) * angle_step;
            }
          }
        }
      }
    };

    struct finalize_alpaka {
      template<typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int chromosome_number,
                                    const int* ligand_num_rotamers,
                                    fp_type* ligand_scores,
                                    fp_type* ligand_best_scores,
                                    chromosome* chromosomes,
                                    chromosome* best_chromosomes) const {
        const int ligand_id = static_cast<int>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        const int thread_id = static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        if (thread_id != 0) {
          return;
        }

        const int num_rotamers = ligand_num_rotamers[ligand_id];
        chromosome* l_chromosomes = chromosomes + ligand_id * chromosome_number;
        fp_type* scores = ligand_scores + chromosome_number * ligand_id;

        int min_index = 0;
        fp_type min_score = scores[0];
        for (int chromosome_index = 1; chromosome_index < chromosome_number; ++chromosome_index) {
          if (min_score > scores[chromosome_index]) {
            min_index = chromosome_index;
            min_score = scores[chromosome_index];
          }
        }

        ligand_best_scores[ligand_id] = min_score;
        fp_type* dst = (best_chromosomes + ligand_id)->data();
        const fp_type* src = (l_chromosomes + min_index)->data();
        for (int i = 0; i < (6 + num_rotamers); ++i) {
          dst[i] = src[i];
        }
      }
    };
  } // namespace

  template<>
  void genetic_kernel<queue_alpaka>::initialize() {
    q->invoke_kernel<initialize_alpaka>(batch_ligands,
                                        MUDOCK_ALPAKA_BLOCK_SIZE,
                                        population_number,
                                        num_rotamers_b,
                                        population,
                                        scores_b,
                                        seed);
  }

  template<>
  void genetic_kernel<queue_alpaka>::operator()() {
    const auto generation_seed = seed++;
    q->invoke_kernel<iterate_alpaka>(batch_ligands,
                                     MUDOCK_ALPAKA_BLOCK_SIZE,
                                     tournament_length,
                                     mutation_prob,
                                     population_number,
                                     num_rotamers_b,
                                     population,
                                     next_population,
                                     scores_b,
                                     generation_seed);
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
