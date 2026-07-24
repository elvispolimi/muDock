#include <mudock/compute/devices_memory.hpp>
#include <mudock/sycl_implementation/genetic_sycl.hpp>
#include <mudock/sycl_implementation/invoke_kernel_sycl.hpp>
#include <mudock/sycl_implementation/sycl_random.hpp>
#include <mudock/utils.hpp>
#include <sycl/sycl.hpp>

#ifndef MUDOCK_SYCL_WG_SIZE
  #define MUDOCK_SYCL_WG_SIZE 32
#endif

namespace mudock {

  device_memory<sycl_random_object>& get_sycl_random_memory() {
    // Intentionally leaked to avoid thread-local destruction after SYCL runtime teardown.
    thread_local auto* storage = new device_memory<sycl_random_object>();
    return *storage;
  }

  static constexpr fp_type coordinate_step = static_cast<fp_type>(0.2);
  static constexpr fp_type angle_step      = static_cast<fp_type>(4);

  // TODO check the real randomness
  template<typename T>
  const T random_gen_sycl(XORWOWState& state, const T min, const T max) {
    fp_type value;
    if constexpr (is_debug())
      // TODO value here for debug
      value = static_cast<fp_type>(0.4);
    else {
      value = state.next();
    }
    return static_cast<T>((value * static_cast<fp_type>(max - min)) + static_cast<fp_type>(min));
  }

  inline int get_selection_distribution(XORWOWState& state, const int* population_number) {
    return random_gen_sycl<int>(state, 0, *population_number - 1);
  };

  inline fp_type get_init_change_distribution(XORWOWState& state) {
    return random_gen_sycl<fp_type>(state, -45, 45);
  }
  inline fp_type get_mutation_change_distribution(XORWOWState& state) {
    return random_gen_sycl<fp_type>(state, -10, 10);
  };
  inline fp_type get_mutation_coin_distribution(XORWOWState& state) {
    return random_gen_sycl<fp_type>(state, 0, 1);
  };
  inline int get_crossover_distribution(XORWOWState& state, const int* num_rotamers) {
    return random_gen_sycl<int>(state, 0, 6 + *num_rotamers);
  };

  inline int tournament_selection_sycl(XORWOWState& state,
                                       const int tournament_length,
                                       const int chromosome_number,
                                       const fp_type* __restrict__ scores) {
    const int num_iterations = tournament_length;
    int best_individual      = get_selection_distribution(state, &chromosome_number);
    for (int i = 0; i < num_iterations; ++i) {
      auto contended = get_selection_distribution(state, &chromosome_number);
      if (scores[contended] < scores[best_individual]) {
        best_individual = contended;
      }
    }
    return best_individual;
  }

  struct initialize_gpu {
    void operator()(sycl::nd_item<3> it,
                    const int chromosome_number,
                    const int* __restrict__ ligand_num_rotamers,
                    chromosome* __restrict__ chromosomes,
                    XORWOWState* __restrict__ state,
                    fp_type* __restrict__ ligand_scores) const {
      const int ligand_id        = static_cast<int>(it.get_group(0));
      const int local_thread_id  = static_cast<int>(it.get_local_id(0));
      const int thread_per_block = static_cast<int>(it.get_local_range(0));
      const int global_thread_id = local_thread_id + thread_per_block * ligand_id;

      const int num_rotamers       = ligand_num_rotamers[ligand_id];
      chromosome* l_chromosomes    = chromosomes + ligand_id * chromosome_number;
      fp_type* __restrict__ scores = ligand_scores + chromosome_number * ligand_id;
      XORWOWState l_state          = (state[global_thread_id]);

      // Shared memory
      // extern __shared__ fp_type shared_data[];
      // fp_type* s_chromosome_scores = shared_data;
      // Initialize shared scores
      for (int chromosome_index = local_thread_id; chromosome_index < chromosome_number;
           chromosome_index += thread_per_block)
        scores[chromosome_index] = std::numeric_limits<fp_type>::infinity(); // Set initial score value

      // Generate initial population
      for (int chromosome_index = local_thread_id; chromosome_index < chromosome_number;
           chromosome_index += thread_per_block) {
        chromosome& chromo = *(l_chromosomes + chromosome_index);
        MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
        for (int i{0}; i < 3; ++i) { // initialize the rigid translation
          chromo[i] = get_init_change_distribution(l_state) * coordinate_step;
        }
        MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
        for (int i{3}; i < 6 + num_rotamers; ++i) { // initialize the rotations
          chromo[i] = get_init_change_distribution(l_state) * angle_step;
        }
      }

      state[global_thread_id] = l_state;
    }
  };

  struct iterate_gpu {
    void operator()(sycl::nd_item<3> it,
                    const int tournament_length,
                    const fp_type mutation_prob,
                    const int chromosome_number,
                    const int* __restrict__ ligand_num_rotamers,
                    chromosome* __restrict__ chromosomes,
                    chromosome* __restrict__ next_chromosomes,
                    XORWOWState* __restrict__ state,
                    fp_type* __restrict__ ligand_scores) const {
      const int ligand_id        = static_cast<int>(it.get_group(0));
      const int local_thread_id  = static_cast<int>(it.get_local_id(0));
      const int thread_per_block = static_cast<int>(it.get_local_range(0));
      const int global_thread_id = local_thread_id + thread_per_block * ligand_id;

      const int num_rotamers                      = ligand_num_rotamers[ligand_id];
      chromosome* __restrict__ l_chromosomes      = chromosomes + ligand_id * chromosome_number;
      chromosome* __restrict__ l_next_chromosomes = next_chromosomes + ligand_id * chromosome_number;
      XORWOWState l_state                         = (state[global_thread_id]);
      const fp_type* __restrict__ scores          = ligand_scores + chromosome_number * ligand_id;

      // Generate the new population
      for (int chromosome_index = local_thread_id; chromosome_index < chromosome_number;
           chromosome_index += thread_per_block) {
        chromosome& next_chromosome = *(l_next_chromosomes + chromosome_index);

        // select the parent
        // TODO check probably they are always the same
        const int best_individual_1 =
            tournament_selection_sycl(l_state, tournament_length, chromosome_number, scores);
        const int best_individual_2 =
            tournament_selection_sycl(l_state, tournament_length, chromosome_number, scores);

        // generate the offspring
        const int split_index = get_crossover_distribution(l_state, &num_rotamers);
        fp_type* dst          = next_chromosome.data();
        const fp_type* p1     = l_chromosomes[best_individual_1].data();
        const fp_type* p2     = l_chromosomes[best_individual_2].data();
        for (int i = 0; i < (6 + num_rotamers); ++i) { dst[i] = (i < split_index) ? p1[i] : p2[i]; }

        // mutate the offspring
        MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
        for (int i{0}; i < 3; ++i) {
          if (get_mutation_coin_distribution(l_state) < mutation_prob)
            next_chromosome[i] += get_mutation_change_distribution(l_state) * coordinate_step;
        }
        MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
        for (int i{3}; i < 6 + num_rotamers; ++i) {
          if (get_mutation_coin_distribution(l_state) < mutation_prob) {
            next_chromosome[i] += get_mutation_change_distribution(l_state) * angle_step;
          }
        }
      }
      state[global_thread_id] = l_state;
    }
  };

  struct finalize_gpu {
    void operator()(sycl::nd_item<3> it,
                    const int chromosome_number,
                    const int* __restrict__ ligand_num_rotamers,
                    fp_type* __restrict__ ligand_scores,
                    fp_type* __restrict__ ligand_best_scores,
                    chromosome* __restrict__ chromosomes,
                    chromosome* __restrict__ best_chromosomes) const {
      const int ligand_id        = static_cast<int>(it.get_group(0));
      const int local_thread_id  = static_cast<int>(it.get_local_id(0));
      const int thread_per_block = static_cast<int>(it.get_local_range(0));
      const auto& sub_group      = it.get_sub_group();

      const int num_rotamers                 = ligand_num_rotamers[ligand_id];
      chromosome* __restrict__ l_chromosomes = chromosomes + ligand_id * chromosome_number;
      fp_type* __restrict__ scores           = ligand_scores + chromosome_number * ligand_id;

      // Compute the maximum value within the warp
      // Assuming each warp has 32 threads
      int min_index = local_thread_id;
      fp_type min_score =
          min_index < chromosome_number ? scores[min_index] : std::numeric_limits<fp_type>::infinity();
      for (int chromosome_index = local_thread_id + thread_per_block; chromosome_index < chromosome_number;
           chromosome_index += thread_per_block) {
        if (min_score > scores[chromosome_index]) {
          min_index = chromosome_index;
          min_score = scores[chromosome_index];
        }
      }
      const fp_type best_score = sycl::reduce_over_group(sub_group, min_score, sycl::minimum());

      // TODO checks that only one can do it
      if (min_score == best_score) {
        ligand_best_scores[ligand_id] = min_score;
        memcpy((*(best_chromosomes + ligand_id)).data(),
               (*(l_chromosomes + min_index)).data(),
               sizeof(fp_type) * (6 + num_rotamers));
      }
    }
  };

  template<>
  void genetic_kernel<queue_sycl>::operator()() {
    auto& random_memory = get_sycl_random_memory();
    q->invoke_kernel<iterate_gpu>(batch_ligands,
                                  MUDOCK_SYCL_WG_SIZE,
                                  tournament_length,
                                  mutation_prob,
                                  population_number,
                                  num_rotamers_b,
                                  population,
                                  next_population,
                                  random_memory.get_data()->dev_pointer(),
                                  scores_b);
  }
  template<>
  void genetic_kernel<queue_sycl>::initialize() {
    auto& random_memory = get_sycl_random_memory();
    // TODO each time or once per computation starts
    random_memory.init(q);
    // TODO chek assumption on num_threads
    random_memory.get_data()->alloc(batch_ligands * MUDOCK_SYCL_WG_SIZE, seed);

    q->invoke_kernel<initialize_gpu>(batch_ligands,
                                     MUDOCK_SYCL_WG_SIZE,
                                     population_number,
                                     num_rotamers_b,
                                     population,
                                     random_memory.get_data()->dev_pointer(),
                                     scores_b);
  }
  template<>
  void genetic_kernel<queue_sycl>::finalize() {
    q->invoke_kernel<finalize_gpu>(batch_ligands,
                                   MUDOCK_SYCL_WG_SIZE,
                                   population_number,
                                   num_rotamers_b,
                                   scores_b,
                                   best_scores_b,
                                   population,
                                   best_chromosomes_b);
  }
} // namespace mudock
