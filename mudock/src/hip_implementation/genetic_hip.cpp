#include <hiprand/hiprand_kernel.h>
#include <mudock/compute/genetic.hpp>
#include <mudock/compute/geometric_transform.hpp>
#include <mudock/compute/pipeline_registry.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/hip_implementation/hip_utils.hpp>
#include <mudock/hip_implementation/mutate.hpp>
#include <mudock/hip_implementation/queue_hip.hpp>
#include <mudock/hip_implementation/scoring_hip.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  static constexpr fp_type coordinate_step{0.2};
  static constexpr fp_type angle_step{4};

  template<typename T>
  __device__ __forceinline__ const T random_gen_hip(hiprandState& state, const T min, const T max) {
    fp_type value;
    if constexpr (is_debug()) {
      // TODO value here for debug
      value = static_cast<fp_type>(0.4);
    } else {
      value = hiprand_uniform(&state);
    }
    return static_cast<T>((value * static_cast<fp_type>(max - min)) + min);
  }

  __device__ __forceinline__ int get_selection_distribution(hiprandState& state,
                                                            const int* population_number) {
    return random_gen_hip<int>(state, 0, *population_number - 1);
  };

  __device__ __forceinline__ fp_type get_init_change_distribution(hiprandState& state) {
    return random_gen_hip<fp_type>(state, -45, 45);
  }
  __device__ __forceinline__ fp_type get_mutation_change_distribution(hiprandState& state) {
    return random_gen_hip<fp_type>(state, -10, 10);
  };
  __device__ __forceinline__ fp_type get_mutation_coin_distribution(hiprandState& state) {
    return random_gen_hip<fp_type>(state, 0, 1);
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
  __device__ __forceinline__ int get_crossover_distribution(hiprandState& state, const int* num_rotamers) {
    return random_gen_hip<int>(state, 0, 6 + *num_rotamers);
  };

  __device__ __forceinline__ int tournament_selection_hip(hiprandState& state,
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
    using scoring_kernel_t                   = typename scoring_backend_t::device_kernel;
    const auto& geometry                     = genetic.geometry;
    const fp_type max_translation_half_width = genetic.max_translation_half_width;
    const int ligand_id                      = blockIdx.x;
    const int lane                           = threadIdx.x;
    const int num_rotamers                   = genetic.num_rotamers[ligand_id];
    const int num_genes                      = 6 + num_rotamers;
    const int state_index                    = ligand_id * BLOCK_SIZE + lane;
    hiprandState local_state;
    hiprand_init(static_cast<unsigned long long>(genetic.seed) + state_index, state_index, 0, &local_state);
    chromosome* population      = genetic.population + ligand_id * genetic.population_size;
    chromosome* next_population = genetic.next_population + ligand_id * genetic.population_size;
    fp_type* scores             = genetic.scores + ligand_id * genetic.population_size;

    for (int index = lane; index < genetic.population_size; index += BLOCK_SIZE) {
      scores[index]     = std::numeric_limits<fp_type>::infinity();
      chromosome& genes = population[index];
      MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
      for (int gene = 0; gene < 3; ++gene)
        genes[gene] = get_init_change_distribution(local_state) * coordinate_step;
      limit_translation<bounded>(genes, max_translation_half_width);
      for (int gene = 3; gene < num_genes; ++gene)
        genes[gene] = get_init_change_distribution(local_state) * angle_step;
    }
    __syncthreads();

    for (int generation = 0; generation < genetic.generations; ++generation) {
      for (int pose = 0; pose < genetic.population_size; ++pose) {
        apply_hip<MAX_ATOMS>(genetic.population_size,
                             geometry.batch_atom_stride,
                             geometry.original_x,
                             geometry.original_y,
                             geometry.original_z,
                             geometry.x_scratch,
                             geometry.y_scratch,
                             geometry.z_scratch,
                             population - ligand_id * genetic.population_size,
                             geometry.fragments,
                             geometry.ligand_fragments_start,
                             geometry.fragments_start_index,
                             geometry.fragments_stop_index,
                             geometry.frag_indices_start,
                             geometry.num_rotamers,
                             geometry.num_atoms,
                             pose);
        __syncthreads();
      }

      for (int pose = 0; pose < genetic.population_size; ++pose) {
        const int coordinate_offset =
            (ligand_id * genetic.population_size + pose) * geometry.batch_atom_stride;
        auto ligand_scoring              = scoring;
        ligand_scoring.ligand_index      = ligand_id;
        ligand_scoring.x                 = geometry.x_scratch + coordinate_offset;
        ligand_scoring.y                 = geometry.y_scratch + coordinate_offset;
        ligand_scoring.z                 = geometry.z_scratch + coordinate_offset;
        ligand_scoring.batch_atom_stride = geometry.batch_atom_stride;
        const fp_type score              = scoring_kernel_t::template score<MAX_ATOMS>(ligand_scoring);
        if (lane == 0)
          scores[pose] = score;
        __syncthreads();
      }

      if (generation + 1 < genetic.generations) {
        for (int index = lane; index < genetic.population_size; index += BLOCK_SIZE) {
          chromosome& child = next_population[index];
          const int parent1 = tournament_selection_hip(local_state,
                                                       genetic.tournament_length,
                                                       genetic.population_size,
                                                       scores);
          const int parent2 = tournament_selection_hip(local_state,
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
          limit_translation<bounded>(child, max_translation_half_width);
          for (int gene = 3; gene < num_genes; ++gene) {
            if (get_mutation_coin_distribution(local_state) < genetic.mutation_probability)
              child[gene] += get_mutation_change_distribution(local_state) * angle_step;
          }
        }
        __syncthreads();
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
      const fp_type other_score = SHFL_DOWN(BITLANE_MASK, best_score, offset, BLOCK_SIZE);
      const int other_index     = SHFL_DOWN(BITLANE_MASK, best_index, offset, BLOCK_SIZE);
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

  template<>
  template<scoring_function scoring>
  void genetic_backend<queue_hip>::run(
      std::shared_ptr<queue_hip> q,
      const genetic_args& genetic,
      geometric<queue_hip>&,
      typename scoring_backend<queue_hip, scoring>::stage_type& scoring_stage) {
    auto scoring_args         = scoring_stage.kernel_args();
    const auto kernel_args    = genetic;
    void* kernel_parameters[] = {(void*) &kernel_args, (void*) &scoring_args};
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          constexpr int max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          using selected_scoring  = scoring_backend<queue_hip, scoring>;
          const auto kernel       = genetic.max_translation_half_width > fp_type{0}
                                        ? (void*) genetic_gpu<max_atoms, true, selected_scoring>
                                        : (void*) genetic_gpu<max_atoms, false, selected_scoring>;
          q->launch_kernel(kernel, kernel_parameters, genetic.batch_ligands, BLOCK_SIZE);
        },
        genetic.geometry.batch_atom_stride,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

  [[maybe_unused]] extern const auto registered_genetic_hip_backends =
      make_genetic_backend_references<genetic_backend<queue_hip>>(genetic_pipeline_registry{});

} // namespace mudock
