#include <mudock/compute/devices_memory.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cuda_implementation/cuda_random.cuh>
#include <mudock/cuda_implementation/genetic_cuda.cuh>
#include <mudock/cuda_implementation/queue_cuda.cuh>
#include <mudock/utils.hpp>

namespace mudock {
  thread_local device_memory<cuda_random_object> cuda_random_memory;

  static constexpr fp_type coordinate_step = static_cast<fp_type>(0.2);
  static constexpr fp_type angle_step{4};

  template<typename T>
  __device__ __forceinline__ const T random_gen_cuda(curandState& state, const T min, const T max) {
    fp_type value;
    if constexpr (is_debug()) {
      // TODO value here for debug
      value = static_cast<fp_type>(0.4);
    } else {
      value = curand_uniform(&state);
    }
    return static_cast<T>((value * static_cast<fp_type>(max - min)) + min);
  }

  __device__ __forceinline__ int get_selection_distribution(curandState& state,
                                                            const int* population_number) {
    return random_gen_cuda<int>(state, 0, *population_number - 1);
  };

  __device__ __forceinline__ fp_type get_init_change_distribution(curandState& state) {
    return random_gen_cuda<fp_type>(state, -45, 45);
  }
  __device__ __forceinline__ fp_type get_mutation_change_distribution(curandState& state) {
    return random_gen_cuda<fp_type>(state, -10, 10);
  };
  __device__ __forceinline__ fp_type get_mutation_coin_distribution(curandState& state) {
    return random_gen_cuda<fp_type>(state, 0, 1);
  };
  // TODO check what happens if max num_rotamers is reached, read for split index could go out of bound
  __device__ __forceinline__ int get_crossover_distribution(curandState& state, const int* num_rotamers) {
    return random_gen_cuda<int>(state, 0, 6 + *num_rotamers);
  };

  __device__ __forceinline__ int tournament_selection_cuda(curandState& state,
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

  __global__ void initialize_gpu(const int chromosome_number,
                                 const int* __restrict__ ligand_num_rotamers,
                                 chromosome* __restrict__ chromosomes,
                                 curandState* __restrict__ state,
                                 fp_type* __restrict__ ligand_scores) {
    const int ligand_id        = blockIdx.x;
    const int local_thread_id  = threadIdx.x;
    const int thread_per_block = blockDim.x;
    const int global_thread_id = local_thread_id + thread_per_block * ligand_id;

    const int num_rotamers       = ligand_num_rotamers[ligand_id];
    chromosome* l_chromosomes    = chromosomes + ligand_id * chromosome_number;
    fp_type* __restrict__ scores = ligand_scores + chromosome_number * ligand_id;
    curandState l_state          = (state[global_thread_id]);

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

  __global__ void iterate_gpu(const int tournament_length,
                              const fp_type mutation_prob,
                              const int chromosome_number,
                              const int* __restrict__ ligand_num_rotamers,
                              chromosome* __restrict__ chromosomes,
                              chromosome* __restrict__ next_chromosomes,
                              curandState* __restrict__ state,
                              fp_type* __restrict__ ligand_scores) {
    const int ligand_id        = blockIdx.x;
    const int local_thread_id  = threadIdx.x;
    const int thread_per_block = blockDim.x;
    const int global_thread_id = local_thread_id + thread_per_block * ligand_id;

    const int num_rotamers                      = ligand_num_rotamers[ligand_id];
    chromosome* __restrict__ l_chromosomes      = chromosomes + ligand_id * chromosome_number;
    chromosome* __restrict__ l_next_chromosomes = next_chromosomes + ligand_id * chromosome_number;
    curandState l_state                         = (state[global_thread_id]);
    const fp_type* __restrict__ scores          = ligand_scores + chromosome_number * ligand_id;

    // Generate the new population
    for (int chromosome_index = local_thread_id; chromosome_index < chromosome_number;
         chromosome_index += thread_per_block) {
      chromosome& next_chromosome = *(l_next_chromosomes + chromosome_index);

      // select the parent
      // TODO check probably they are always the same
      const int best_individual_1 =
          tournament_selection_cuda(l_state, tournament_length, chromosome_number, scores);
      const int best_individual_2 =
          tournament_selection_cuda(l_state, tournament_length, chromosome_number, scores);

      // generate the offspring
      const int split_index = get_crossover_distribution(l_state, &num_rotamers);
      // memcpy(next_chromosome.data(), &(l_chromosomes[best_individual_1][0]), split_index * sizeof(fp_type));
      // const int parent2_copy_size = 6 + num_rotamers - split_index;
      // if (parent2_copy_size > 0)
      //   memcpy(next_chromosome.data() + split_index,
      //          &(l_chromosomes[best_individual_2][split_index]),
      //          parent2_copy_size * sizeof(fp_type));
      fp_type* dst      = next_chromosome.data();
      const fp_type* p1 = l_chromosomes[best_individual_1].data();
      const fp_type* p2 = l_chromosomes[best_individual_2].data();
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

  __global__ void finalize_gpu(const int chromosome_number,
                               const int* __restrict__ ligand_num_rotamers,
                               fp_type* __restrict__ ligand_scores,
                               fp_type* __restrict__ ligand_best_scores,
                               chromosome* __restrict__ chromosomes,
                               chromosome* __restrict__ best_chromosomes) {
    const int ligand_id        = blockIdx.x;
    const int local_thread_id  = threadIdx.x;
    const int thread_per_block = blockDim.x;

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
    // Intra warp reduction
    MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
    for (int offset = BLOCK_SIZE / 2; offset > 0; offset /= 2) {
      const fp_type other_min_score = __shfl_down_sync(0xFFFFFFFF, min_score, offset);
      const int other_min_index     = __shfl_down_sync(0xFFFFFFFF, min_index, offset);
      if (other_min_score < min_score) {
        min_score = other_min_score;
        min_index = other_min_index;
      }
    }
    if (local_thread_id == 0) {
      ligand_best_scores[ligand_id] = min_score;

      fp_type* dst       = (best_chromosomes + ligand_id)->data();
      const fp_type* src = (l_chromosomes + min_index)->data();
      for (int i = 0; i < (6 + num_rotamers); ++i) { dst[i] = src[i]; }
    }
  }

  template<>
  void genetic_kernel<queue_cuda>::initialize() {
    // TODO each time or once per computation starts
    cuda_random_memory.init(q);
    // TODO chek assumption on num_threads
    cuda_random_memory.get_data()->alloc(batch_ligands * BLOCK_SIZE, seed);

    void* args[] = {(void*) &population_number,
                    (void*) &num_rotamers_b,
                    (void*) &population,
                    (void*) cuda_random_memory.get_data()->dev_pointer_ref(),
                    (void*) &scores_b};
    //TODO check grid/dimensions
    q->launch_kernel((void*) initialize_gpu, args, batch_ligands, BLOCK_SIZE);
  }
  template<>
  void genetic_kernel<queue_cuda>::operator()() {
    void* args[] = {(void*) &tournament_length,
                    (void*) &mutation_prob,
                    (void*) &population_number,
                    (void*) &num_rotamers_b,
                    (void*) &population,
                    (void*) &next_population,
                    (void*) cuda_random_memory.get_data()->dev_pointer_ref(),
                    (void*) &scores_b};
    //TODO check grid/dimensions
    q->launch_kernel((void*) iterate_gpu, args, batch_ligands, BLOCK_SIZE);
  }
  template<>
  void genetic_kernel<queue_cuda>::finalize() {
    void* args[] = {(void*) &population_number,
                    (void*) &num_rotamers_b,
                    (void*) &scores_b,
                    (void*) &best_scores_b,
                    (void*) &population,
                    (void*) &best_chromosomes_b};
    //TODO check grid/dimensions
    q->launch_kernel((void*) finalize_gpu, args, batch_ligands, BLOCK_SIZE);
  }
} // namespace mudock
