#pragma once

#include <hiprand/hiprand.h>
#include <hiprand/hiprand_kernel.h>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/grid.hpp>
#include <mudock/hip_implementation/calc_energy.hpp>
#include <mudock/hip_implementation/hip_utils.hpp>
#include <mudock/hip_implementation/mutate.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  static constexpr fp_type coordinate_step{0.2};
  static constexpr fp_type angle_step{4};

  template<typename T>
  __device__ inline const T random_gen_hip(hiprandState& state, const T min, const T max) {
    fp_type value;
    if constexpr (is_debug()) {
      // TODO value here for debug
      value = fp_type{0.4};
    } else {
      value = hiprand_uniform(&state);
    }
    return static_cast<T>((value * static_cast<fp_type>(max - min)) + min);
  }

  __device__ inline int get_selection_distribution(hiprandState& state, const int* population_number) {
    return random_gen_hip<int>(state, 0, *population_number - 1);
  };

  __device__ inline fp_type get_init_change_distribution(hiprandState& state) {
    return random_gen_hip<fp_type>(state, -45, 45);
  }
  __device__ inline fp_type get_mutation_change_distribution(hiprandState& state) {
    return random_gen_hip<fp_type>(state, -10, 10);
  };
  __device__ inline fp_type get_mutation_coin_distribution(hiprandState& state) {
    return random_gen_hip<fp_type>(state, 0, 1);
  };
  __device__ inline int get_crossover_distribution(hiprandState& state, const int* num_rotamers) {
    return random_gen_hip<int>(state, 0, 6 + *num_rotamers);
  };

  __device__ inline int tournament_selection_hip(hiprandState& state,
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

  // TODO check the syncwarp
  // TODO OPT: template parameter based on number of atoms, rotamers, chromosomes and population
  // TODO missing bucketizer -> see CUDA version
  // TODO check 0xFFFFFFFF on AMD GPUs
  template<int MAX_ATOMS>
  __global__ void evaluate_fitness(const int num_generations,
                                   const int tournament_length,
                                   const fp_type mutation_prob,
                                   const int chromosome_number,
                                   const int chromosome_stride,
                                   const int atom_stride,
                                   const int map_index_x,
                                   const int map_index_xy,
                                   const int map_index_xyz,
                                   const fp_type* __restrict__ original_ligand_x,
                                   const fp_type* __restrict__ original_ligand_y,
                                   const fp_type* __restrict__ original_ligand_z,
                                   fp_type* __restrict__ scratch_ligand_x,
                                   fp_type* __restrict__ scratch_ligand_y,
                                   fp_type* __restrict__ scratch_ligand_z,
                                   const fp_type* __restrict__ ligand_vol,
                                   const fp_type* __restrict__ ligand_solpar,
                                   const fp_type* __restrict__ ligand_charge,
                                   const int* __restrict__ ligand_num_nonbonds,
                                   const int* __restrict__ ligand_nonbond_a1,
                                   const int* __restrict__ ligand_nonbond_a2,
                                   const fp_type* __restrict__ ligand_nonbond_cA,
                                   const fp_type* __restrict__ ligand_nonbond_cB,
                                   const int* __restrict__ ligand_nonbond_xB,
                                   const int* __restrict__ ligand_num_atoms,
                                   const int* __restrict__ ligand_num_rotamers,
                                   const int* __restrict__ ligand_fragments,
                                   const int* __restrict__ ligand_fragments_start,
                                   const int* __restrict__ frag_start_atom_index,
                                   const int* __restrict__ frag_stop_atom_index,
                                   const int* __restrict__ frag_indices_start,
                                   chromosome* __restrict__ chromosomes,
                                   const fp_type* __restrict__ grid_maps,
                                   const int* __restrict__ map_ligand_offsets,
                                   hiprandState* __restrict__ state,
                                   fp_type* __restrict__ ligand_scores,
                                   chromosome* __restrict__ best_chromosomes) {
    const int ligand_id        = blockIdx.x;
    const int local_thread_id  = threadIdx.x;
    const int thread_per_block = blockDim.x;
    const int global_thread_id = local_thread_id + thread_per_block * ligand_id;

    const int num_atoms    = ligand_num_atoms[ligand_id];
    const int num_nonbonds = ligand_num_nonbonds[ligand_id + 1] - ligand_num_nonbonds[ligand_id];
    const int num_rotamers = ligand_num_rotamers[ligand_id];

    const fp_type* l_original_ligand_x = original_ligand_x + ligand_id * atom_stride;
    const fp_type* l_original_ligand_y = original_ligand_y + ligand_id * atom_stride;
    const fp_type* l_original_ligand_z = original_ligand_z + ligand_id * atom_stride;
    fp_type* l_scratch_ligand_x        = scratch_ligand_x + ligand_id * atom_stride;
    fp_type* l_scratch_ligand_y        = scratch_ligand_y + ligand_id * atom_stride;
    fp_type* l_scratch_ligand_z        = scratch_ligand_z + ligand_id * atom_stride;
    const fp_type* l_ligand_vol        = ligand_vol + ligand_id * atom_stride;
    const fp_type* l_ligand_solpar     = ligand_solpar + ligand_id * atom_stride;
    const fp_type* l_ligand_charge     = ligand_charge + ligand_id * atom_stride;
    chromosome* l_chromosomes          = chromosomes + ligand_id * chromosome_stride;
    // Point to the next population buffer
    chromosome* l_next_chromosomes      = chromosomes + ligand_id * chromosome_stride + chromosome_number;
    const auto* l_fragments             = ligand_fragments + ligand_fragments_start[ligand_id];
    const auto* l_frag_start_atom_index = frag_start_atom_index + frag_indices_start[ligand_id];
    const auto* l_frag_stop_atom_index  = frag_stop_atom_index + frag_indices_start[ligand_id];
    const auto* l_map_ligand_offsets    = map_ligand_offsets + ligand_id * atom_stride;
    const int* l_ligand_nonbond_a1      = ligand_nonbond_a1 + ligand_num_nonbonds[ligand_id];
    const int* l_ligand_nonbond_a2      = ligand_nonbond_a2 + ligand_num_nonbonds[ligand_id];
    const fp_type* l_ligand_nonbond_cA  = ligand_nonbond_cA + ligand_num_nonbonds[ligand_id];
    const fp_type* l_ligand_nonbond_cB  = ligand_nonbond_cB + ligand_num_nonbonds[ligand_id];
    const int* l_ligand_nonbond_xB      = ligand_nonbond_xB + ligand_num_nonbonds[ligand_id];

    hiprandState& l_state = (state[global_thread_id]);

    // Shared memory
    extern __shared__ fp_type shared_data[];
    fp_type* s_chromosome_scores = shared_data;
    // Initialize shared scores
    for (int chromosome_index = local_thread_id; chromosome_index < max(chromosome_number, thread_per_block);
         chromosome_index += thread_per_block)
      s_chromosome_scores[chromosome_index] =
          std::numeric_limits<fp_type>::infinity(); // Set initial score value

    // Generate initial population
    for (int chromosome_index = local_thread_id; chromosome_index < chromosome_number;
         chromosome_index += thread_per_block) {
      chromosome& chromo = *(l_chromosomes + chromosome_index);
#pragma unroll
      for (int i{0}; i < 3; ++i) { // initialize the rigid translation
        chromo[i] = get_init_change_distribution(l_state) * coordinate_step;
      }
#pragma unroll
      for (int i{3}; i < 6 + num_rotamers; ++i) { // initialize the rotations
        chromo[i] = get_init_change_distribution(l_state) * angle_step;
      }
    }
    __syncthreads();
    // TODO maybe template parameter?
    for (int generation = 0; generation < num_generations; ++generation) {
      for (int chromosome_index = 0; chromosome_index < chromosome_number; ++chromosome_index) {
        // Copy original coordinates
        // TODO OPT: shared memory for coordinate ?
        for (int atom_index = local_thread_id; atom_index < MAX_ATOMS; atom_index += thread_per_block) {
          if (atom_index < num_atoms) {
            l_scratch_ligand_x[atom_index] = l_original_ligand_x[atom_index];
            l_scratch_ligand_y[atom_index] = l_original_ligand_y[atom_index];
            l_scratch_ligand_z[atom_index] = l_original_ligand_z[atom_index];
          }
        }
        // Modify coordinates
        apply_hip<MAX_ATOMS>(l_scratch_ligand_x,
                             l_scratch_ligand_y,
                             l_scratch_ligand_z,
                             *(l_chromosomes + chromosome_index),
                             l_fragments,
                             l_frag_start_atom_index,
                             l_frag_stop_atom_index,
                             num_rotamers,
                             num_atoms);

        fp_type energy = calc_energy<MAX_ATOMS>(l_scratch_ligand_x,
                                                l_scratch_ligand_y,
                                                l_scratch_ligand_z,
                                                l_ligand_vol,
                                                l_ligand_solpar,
                                                l_ligand_charge,
                                                num_atoms,
                                                num_rotamers,
                                                num_nonbonds,
                                                l_ligand_nonbond_a1,
                                                l_ligand_nonbond_a2,
                                                l_ligand_nonbond_cA,
                                                l_ligand_nonbond_cB,
                                                l_ligand_nonbond_xB,
                                                map_index_x,
                                                map_index_xy,
                                                map_index_xyz,
                                                grid_maps,
                                                l_map_ligand_offsets);

        // Perform a tree reduction using __shfl_down_sync
        // TODO check performance
        for (int offset = warpSize / 2; offset > 0; offset /= 2) {
          energy += __shfl_down_sync(BITLANE_MASK, energy, offset);
        }

        if (local_thread_id == 0) {
          s_chromosome_scores[chromosome_index] = energy;
        }
      }

      // Generate the new population
      for (int chromosome_index = local_thread_id; chromosome_index < chromosome_number;
           chromosome_index += thread_per_block) {
        chromosome& next_chromosome = *(l_next_chromosomes + chromosome_index);

        // select the parent
        const int best_individual_1 =
            tournament_selection_hip(l_state, tournament_length, chromosome_number, s_chromosome_scores);
        const int best_individual_2 =
            tournament_selection_hip(l_state, tournament_length, chromosome_number, s_chromosome_scores);

        // generate the offspring
        const int split_index = get_crossover_distribution(l_state, &num_rotamers);
        memcpy(next_chromosome.data(), &(l_chromosomes[best_individual_1][0]), split_index * sizeof(fp_type));
        const int parent2_copy_size = 6 + num_rotamers - split_index;
        if (parent2_copy_size > 0)
          memcpy(next_chromosome.data() + split_index,
                 &(l_chromosomes[best_individual_2][split_index]),
                 parent2_copy_size * sizeof(fp_type));

// mutate the offspring
#pragma unroll
        for (int i{0}; i < 3; ++i) {
          if (get_mutation_coin_distribution(l_state) < mutation_prob)
            next_chromosome[i] += get_mutation_change_distribution(l_state) * coordinate_step;
        }
#pragma unroll
        for (int i{3}; i < 6 + num_rotamers; ++i) {
          if (get_mutation_coin_distribution(l_state) < mutation_prob) {
            next_chromosome[i] += get_mutation_change_distribution(l_state) * angle_step;
          }
        }
      }

      // Swap the actual population with the next one
      // TODO check const
      chromosome* const tmp_chromosomes = l_chromosomes;
      l_chromosomes                     = l_next_chromosomes;
      l_next_chromosomes                = tmp_chromosomes;
      __syncthreads();
    }

    // Output

    // Compute the maximum value within the warp
    // Assuming each warp has 32 threads
    int min_index     = local_thread_id;
    fp_type min_score = s_chromosome_scores[min_index];
    for (int chromosome_index = local_thread_id + thread_per_block; chromosome_index < chromosome_number;
         chromosome_index += thread_per_block) {
      if (min_score > s_chromosome_scores[chromosome_index]) {
        min_index = chromosome_index;
        min_score = s_chromosome_scores[chromosome_index];
      }
    }
    // Intra warp reduction
    for (int offset = warpSize / 2; offset > 0; offset /= 2) {
      fp_type other_min_score = __shfl_down_sync(BITLANE_MASK, min_score, offset);
      int other_min_index     = __shfl_down_sync(BITLANE_MASK, min_index, offset);
      if (other_min_score < min_score) {
        min_score = other_min_score;
        min_index = other_min_index;
      }
    }
    if (local_thread_id == 0) {
      const fp_type tors_free_energy = num_rotamers * autodock_parameters::coeff_tors;
      ligand_scores[ligand_id]       = min_score + tors_free_energy;
      memcpy((*(best_chromosomes + ligand_id)).data(),
             (*(l_chromosomes + min_index)).data(),
             sizeof(fp_type) * (6 + num_rotamers));
    }
  }
} // namespace mudock
