#include <curand_kernel.h>
#include <mudock/compute/genetic.hpp>
#include <mudock/compute/genetic_args.hpp>
#include <mudock/compute/pipeline_registry.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>
#include <mudock/cuda_implementation/mutate.cuh>
#include <mudock/cuda_implementation/scoring_cuda.cuh>
#include <mudock/utils.hpp>

namespace mudock {
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

  template<bool bounded>
  __device__ __forceinline__ void limit_translation(chromosome& chrom, const fp_type half_width) {
    if constexpr (!bounded)
      return;
    chrom[0] = fminf(fmaxf(chrom[0], -half_width), half_width);
    chrom[1] = fminf(fmaxf(chrom[1], -half_width), half_width);
    chrom[2] = fminf(fmaxf(chrom[2], -half_width), half_width);
  }
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

  template<int MAX_ATOMS, bool bounded, class scoring_backend_t>
  __global__ void genetic_gpu(genetic_args genetic, typename scoring_backend_t::arguments scoring) {
    using scoring_kernel_t              = typename scoring_backend_t::device_kernel;
    const geom_transform_args& geometry = genetic.geometry;
    const int ligand_id                 = blockIdx.x;
    const int lane                      = threadIdx.x;
    const int num_rotamers              = genetic.num_rotamers[ligand_id];
    const int num_genes                 = 6 + num_rotamers;
    const int state_index               = ligand_id * BLOCK_SIZE + lane;
    curandState local_state;
    curand_init(static_cast<unsigned long long>(genetic.seed) + state_index, state_index, 0, &local_state);
    chromosome* population      = genetic.population + ligand_id * genetic.population_size;
    chromosome* next_population = genetic.next_population + ligand_id * genetic.population_size;
    fp_type* scores             = genetic.scores + ligand_id * genetic.population_size;

    for (int index = lane; index < genetic.population_size; index += BLOCK_SIZE) {
      scores[index]     = std::numeric_limits<fp_type>::infinity();
      chromosome& genes = population[index];
      MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
      for (int gene = 0; gene < 3; ++gene)
        genes[gene] = get_init_change_distribution(local_state) * coordinate_step;
      limit_translation<bounded>(genes, genetic.max_translation_half_width);
      for (int gene = 3; gene < num_genes; ++gene)
        genes[gene] = get_init_change_distribution(local_state) * angle_step;
    }
    __syncwarp();

    for (int generation = 0; generation < genetic.generations; ++generation) {
      for (int pose = 0; pose < genetic.population_size; ++pose) {
        const int coordinate_offset =
            (ligand_id * genetic.population_size + pose) * geometry.batch_atom_stride;
        apply_chromosome_cuda<MAX_ATOMS>(geometry,
                                         ligand_id,
                                         population[pose],
                                         geometry.x_scratch + coordinate_offset,
                                         geometry.y_scratch + coordinate_offset,
                                         geometry.z_scratch + coordinate_offset);
        __syncwarp();
      }

      for (int pose = 0; pose < genetic.population_size; ++pose) {
        const int coordinate_offset =
            (ligand_id * genetic.population_size + pose) * geometry.batch_atom_stride;
        auto ligand_scoring         = scoring;
        ligand_scoring.ligand_index = ligand_id;
        ligand_scoring.x            = scoring.x + coordinate_offset;
        ligand_scoring.y            = scoring.y + coordinate_offset;
        ligand_scoring.z            = scoring.z + coordinate_offset;
        const fp_type score         = scoring_kernel_t::template score<MAX_ATOMS>(ligand_scoring);
        if (lane == 0)
          scores[pose] = score;
        __syncwarp();
      }

      if (generation + 1 < genetic.generations) {
        for (int index = lane; index < genetic.population_size; index += BLOCK_SIZE) {
          chromosome& child = next_population[index];
          const int parent1 = tournament_selection_cuda(local_state,
                                                        genetic.tournament_length,
                                                        genetic.population_size,
                                                        scores);
          const int parent2 = tournament_selection_cuda(local_state,
                                                        genetic.tournament_length,
                                                        genetic.population_size,
                                                        scores);
          const int split   = get_crossover_distribution(local_state, &num_rotamers);
          const fp_type* p1 = population[parent1].data();
          const fp_type* p2 = population[parent2].data();
          for (int gene = 0; gene < num_genes; ++gene) child[gene] = gene < split ? p1[gene] : p2[gene];
          MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
          for (int gene = 0; gene < 3; ++gene) {
            if (get_mutation_coin_distribution(local_state) < genetic.mutation_probability)
              child[gene] += get_mutation_change_distribution(local_state) * coordinate_step;
          }
          limit_translation<bounded>(child, genetic.max_translation_half_width);
          for (int gene = 3; gene < num_genes; ++gene) {
            if (get_mutation_coin_distribution(local_state) < genetic.mutation_probability)
              child[gene] += get_mutation_change_distribution(local_state) * angle_step;
          }
        }
        __syncwarp();
        chromosome* swap = population;
        population       = next_population;
        next_population  = swap;
      }
    }
    int best_index = lane;
    fp_type best_score =
        best_index < genetic.population_size ? scores[best_index] : std::numeric_limits<fp_type>::infinity();
    for (int index = lane + BLOCK_SIZE; index < genetic.population_size; index += BLOCK_SIZE) {
      if (scores[index] < best_score) {
        best_score = scores[index];
        best_index = index;
      }
    }
    for (int offset = BLOCK_SIZE / 2; offset > 0; offset /= 2) {
      const fp_type other_score = __shfl_down_sync(0xffffffff, best_score, offset, BLOCK_SIZE);
      const int other_index     = __shfl_down_sync(0xffffffff, best_index, offset, BLOCK_SIZE);
      if (other_score < best_score) {
        best_score = other_score;
        best_index = other_index;
      }
    }
    if (lane == 0) {
      genetic.best_scores[ligand_id] = best_score;
      fp_type* target                = genetic.best_chromosomes[ligand_id].data();
      const fp_type* source          = population[best_index].data();
      for (int gene = 0; gene < num_genes; ++gene) target[gene] = source[gene];
    }
  }

  template<typename queue_t>
  template<scoring_function scoring>
  void genetic_backend<queue_t>::run(std::shared_ptr<queue_t> q,
                                     const genetic_args& genetic,
                                     geometric<queue_t>&,
                                     typename scoring_backend<queue_t, scoring>::stage_type& scoring_stage) {
    const auto kernel_args    = genetic;
    const auto scoring_args   = scoring_stage.kernel_args();
    void* kernel_parameters[] = {(void*) &kernel_args, (void*) &scoring_args};

    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          constexpr int max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          using selected_scoring  = scoring_backend<queue_t, scoring>;
          const auto kernel       = kernel_args.max_translation_half_width > fp_type{0}
                                        ? (void*) genetic_gpu<max_atoms, true, selected_scoring>
                                        : (void*) genetic_gpu<max_atoms, false, selected_scoring>;
          q->launch_kernel(kernel, kernel_parameters, kernel_args.batch_ligands, BLOCK_SIZE);
        },
        kernel_args.geometry.batch_atom_stride,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

  [[maybe_unused]] extern const auto registered_genetic_backends =
      make_genetic_backend_references<genetic_backend<queue_cuda>>(genetic_pipeline_registry{});
} // namespace mudock
