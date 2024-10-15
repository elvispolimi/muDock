#include <cstdio>
#include <curand_kernel.h>
#include <mudock/cuda_implementation/calc_energy.cuh>
#include <mudock/cuda_implementation/cuda_check_error_macro.cuh>
#include <mudock/cuda_implementation/evaluate_fitness.cuh>
#include <mudock/cuda_implementation/mutate.cuh>
#include <mudock/grid.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  static constexpr fp_type coordinate_step{0.2};
  static constexpr fp_type angle_step{4};

  __device__ __constant__ fp_type map_min_const[3];
  __device__ __constant__ fp_type map_max_const[3];
  __device__ __constant__ fp_type map_center_const[3];

  void setup_constant_memory(const point3D& minimum_coord,
                             const point3D& maximum_coord,
                             const point3D& center) {
    const fp_type l_map_min[3]{minimum_coord.x, minimum_coord.y, minimum_coord.z};
    const fp_type l_map_max[3]{maximum_coord.x, maximum_coord.y, maximum_coord.z};
    const fp_type l_map_center[3]{center.x, center.y, center.z};

    MUDOCK_CHECK(
        cudaMemcpyToSymbol(map_min_const, &l_map_min, 3 * sizeof(fp_type), 0, cudaMemcpyHostToDevice));
    MUDOCK_CHECK(
        cudaMemcpyToSymbol(map_max_const, &l_map_max, 3 * sizeof(fp_type), 0, cudaMemcpyHostToDevice));
    MUDOCK_CHECK(
        cudaMemcpyToSymbol(map_center_const, &l_map_center, 3 * sizeof(fp_type), 0, cudaMemcpyHostToDevice));
  }

  __device__ fp_type trilinear_interpolation_cuda(const fp_type coord[], const cudaTextureObject_t& tex) {
    // Interpolation CUDA
    const int u0      = coord[0];
    const fp_type p0u = coord[0] - static_cast<fp_type>(u0);
    const fp_type p1u = fp_type{1} - p0u;

    const int v0      = coord[1];
    const fp_type p0v = coord[1] - static_cast<fp_type>(v0);
    const fp_type p1v = fp_type{1} - p0v;

    const int w0      = coord[2];
    const fp_type p0w = coord[2] - static_cast<fp_type>(w0);
    const fp_type p1w = fp_type{1} - p0w;

    const fp_type pu[2] = {p1u, p0u};
    const fp_type pv[2] = {p1v, p0v};
    const fp_type pw[2] = {p1w, p0w};
    fp_type value{0};
#pragma unroll
    for (int i = 0; i <= 1; i++)
#pragma unroll
      for (int t = 0; t <= 1; t++)
#pragma unroll
        for (int n = 0; n <= 1; n++) {
          const fp_type tmp = tex3D<fp_type>(tex, u0 + n, v0 + t, w0 + i);
          value += pu[n] * pv[t] * pw[i] * tmp;
        }
    return value;
  }

  template<typename T>
  __device__ const T random_gen_cuda(curandState& state, const T min, const T max) {
    fp_type value;
    if constexpr (is_debug()) {
      // TODO value here for debug
      value = fp_type{0.4};
    } else {
      value = curand_uniform(&state);
    }
    return static_cast<T>((value * static_cast<fp_type>(max - min)) + min);
  }

  __device__ int get_selection_distribution(curandState& state, const int* population_number) {
    return random_gen_cuda<int>(state, 0, *population_number - 1);
  };

  __device__ fp_type get_init_change_distribution(curandState& state) {
    return random_gen_cuda<fp_type>(state, -45, 45);
  }
  __device__ fp_type get_mutation_change_distribution(curandState& state) {
    return random_gen_cuda<fp_type>(state, -10, 10);
  };
  __device__ fp_type get_mutation_coin_distribution(curandState& state) {
    return random_gen_cuda<fp_type>(state, 0, 1);
  };
  __device__ int get_crossover_distribution(curandState& state, const int* num_rotamers) {
    return random_gen_cuda<int>(state, 0, 6 + *num_rotamers);
  };

  __device__ int tournament_selection_cuda(curandState& state,
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
  //TODO OPT: template parameter based on number of atoms, rotamers, chromosomes and population
  // Interesting the usage of the bucketizer
  __global__ void evaluate_fitness(const int num_generations,
                                   const int tournament_length,
                                   const fp_type mutation_prob,
                                   const int chromosome_number,
                                   const int chromosome_stride,
                                   const int atom_stride,
                                   const int rotamers_stride,
                                   const int nonbond_stride,
                                   const fp_type* __restrict__ original_ligand_x,
                                   const fp_type* __restrict__ original_ligand_y,
                                   const fp_type* __restrict__ original_ligand_z,
                                   fp_type* __restrict__ scratch_ligand_x,
                                   fp_type* __restrict__ scratch_ligand_y,
                                   fp_type* __restrict__ scratch_ligand_z,
                                   const fp_type* __restrict__ ligand_vol,
                                   const fp_type* __restrict__ ligand_solpar,
                                   const fp_type* __restrict__ ligand_charge,
                                   const int* __restrict__ ligand_num_hbond,
                                   const fp_type* __restrict__ ligand_Rij_hb,
                                   const fp_type* __restrict__ ligand_Rii,
                                   const fp_type* __restrict__ ligand_epsij_hb,
                                   const fp_type* __restrict__ ligand_epsii,
                                   const int* __restrict__ ligand_num_nonbonds,
                                   const int* __restrict__ ligand_nonbond_a1,
                                   const int* __restrict__ ligand_nonbond_a2,
                                   const int* __restrict__ ligand_num_atoms,
                                   const int* __restrict__ ligand_num_rotamers,
                                   const int* __restrict__ ligand_fragments,
                                   const int* __restrict__ frag_start_atom_index,
                                   const int* __restrict__ frag_stop_atom_index,
                                   chromosome* __restrict__ chromosomes,
                                   const cudaTextureObject_t* __restrict__ atom_textures,
                                   const int* __restrict__ atom_tex_indexes,
                                   const cudaTextureObject_t electro_texture,
                                   const cudaTextureObject_t desolv_texture,
                                   curandState* __restrict__ state,
                                   fp_type* __restrict__ ligand_scores,
                                   chromosome* __restrict__ best_chromosomes) {
    const int ligand_id        = blockIdx.x;
    const int local_thread_id  = threadIdx.x;
    const int thread_per_block = blockDim.x;
    const int global_thread_id = local_thread_id + thread_per_block * ligand_id;

    const int num_atoms    = ligand_num_atoms[ligand_id];
    const int num_nonbonds = ligand_num_nonbonds[ligand_id];
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
    const int* l_ligand_num_hbond      = ligand_num_hbond + ligand_id * atom_stride;
    const fp_type* l_ligand_Rij_hb     = ligand_Rij_hb + ligand_id * atom_stride;
    const fp_type* l_ligand_Rii        = ligand_Rii + ligand_id * atom_stride;
    const fp_type* l_ligand_epsij_hb   = ligand_epsij_hb + ligand_id * atom_stride;
    const fp_type* l_ligand_epsii      = ligand_epsii + ligand_id * atom_stride;
    chromosome* l_chromosomes          = chromosomes + ligand_id * chromosome_stride;
    // Point to the next population buffer
    chromosome* l_next_chromosomes      = chromosomes + ligand_id * chromosome_stride + chromosome_number;
    const auto* l_fragments             = ligand_fragments + ligand_id * atom_stride * rotamers_stride;
    const auto* l_frag_start_atom_index = frag_start_atom_index + ligand_id * rotamers_stride;
    const auto* l_frag_stop_atom_index  = frag_stop_atom_index + ligand_id * rotamers_stride;
    const auto* l_atom_tex_indexes      = atom_tex_indexes + ligand_id * atom_stride;
    const int* l_ligand_nonbond_a1      = ligand_nonbond_a1 + ligand_id * nonbond_stride;
    const int* l_ligand_nonbond_a2      = ligand_nonbond_a2 + ligand_id * nonbond_stride;
    curandState& l_state                = (state[global_thread_id]);

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
    __syncwarp();
    // TODO maybe template parameter?
    for (int generation = 0; generation < num_generations; ++generation) {
      for (int chromosome_index = 0; chromosome_index < chromosome_number; ++chromosome_index) {
        // Copy original coordinates
        // TODO OPT: shared memory for coordinate ?
        for (int atom_index = local_thread_id; atom_index < num_atoms; atom_index += thread_per_block) {
          l_scratch_ligand_x[atom_index] = l_original_ligand_x[atom_index];
          l_scratch_ligand_y[atom_index] = l_original_ligand_y[atom_index];
          l_scratch_ligand_z[atom_index] = l_original_ligand_z[atom_index];
        }
        // Modify coordinates
        apply_cuda(l_scratch_ligand_x,
                   l_scratch_ligand_y,
                   l_scratch_ligand_z,
                   *(l_chromosomes + chromosome_index),
                   l_fragments,
                   l_frag_start_atom_index,
                   l_frag_stop_atom_index,
                   num_rotamers,
                   atom_stride,
                   num_atoms);

        // Calculate energy
        fp_type elect_total_trilinear = 0;
        fp_type emap_total_trilinear  = 0;
        fp_type dmap_total_trilinear  = 0;
        for (int atom_index = local_thread_id; atom_index < num_atoms; atom_index += thread_per_block) {
          fp_type coord_tex[3]{l_scratch_ligand_x[atom_index],
                               l_scratch_ligand_y[atom_index],
                               l_scratch_ligand_z[atom_index]};

          if (coord_tex[0] < map_min_const[0] || coord_tex[0] > map_max_const[0] ||
              coord_tex[1] < map_min_const[1] || coord_tex[1] > map_max_const[1] ||
              coord_tex[2] < map_min_const[2] || coord_tex[2] > map_max_const[2]) {
            // Is outside
            const fp_type distance_two = powf(fabs(coord_tex[0] - map_center_const[0]), fp_type{2}) +
                                         powf(fabs(coord_tex[1] - map_center_const[1]), fp_type{2}) +
                                         powf(fabs(coord_tex[2] - map_center_const[2]), fp_type{2});

            const fp_type epenalty = distance_two * ENERGYPENALTY;
            elect_total_trilinear += epenalty;
            emap_total_trilinear += epenalty;
          } else {
            // Is inside
            // Center atom coordinates on the grid center
            coord_tex[0] = (l_scratch_ligand_x[atom_index] - map_min_const[0]) * inv_spacing,
            coord_tex[1] = (l_scratch_ligand_y[atom_index] - map_min_const[1]) * inv_spacing;
            coord_tex[2] = (l_scratch_ligand_z[atom_index] - map_min_const[2]) * inv_spacing;
            //  TODO check approximations with in hardware interpolation
            // elect_total_trilinear +=
            //     tex3D<fp_type>(electro_texture, coord_tex[0], coord_tex[1], coord_tex[2]) *
            //     l_ligand_charge[atom_index];
            // dmap_total_trilinear += tex3D<fp_type>(desolv_texture, coord_tex[0], coord_tex[1], coord_tex[2]) *
            //                         fabsf(l_ligand_charge[atom_index]);
            // emap_total_trilinear += tex3D<fp_type>(atom_textures[l_atom_tex_indexes[atom_index]],
            //                                        coord_tex[0],
            //                                        coord_tex[1],
            //                                        coord_tex[2]);

            elect_total_trilinear +=
                trilinear_interpolation_cuda(coord_tex, electro_texture) * l_ligand_charge[atom_index];
            dmap_total_trilinear +=
                trilinear_interpolation_cuda(coord_tex, desolv_texture) * fabsf(l_ligand_charge[atom_index]);
            emap_total_trilinear +=
                trilinear_interpolation_cuda(coord_tex, atom_textures[l_atom_tex_indexes[atom_index]]);
          }
        }
        fp_type total_trilinear = elect_total_trilinear + dmap_total_trilinear + emap_total_trilinear;

        fp_type total_eintcal{0};
        if (num_rotamers > 0)
          total_eintcal += calc_intra_energy(l_scratch_ligand_x,
                                             l_scratch_ligand_y,
                                             l_scratch_ligand_z,
                                             l_ligand_vol,
                                             l_ligand_solpar,
                                             l_ligand_charge,
                                             l_ligand_num_hbond,
                                             l_ligand_Rij_hb,
                                             l_ligand_Rii,
                                             l_ligand_epsij_hb,
                                             l_ligand_epsii,
                                             num_nonbonds,
                                             l_ligand_nonbond_a1,
                                             l_ligand_nonbond_a2);

        // Perform a tree reduction using __shfl_down_sync
        // TODO check performance
        for (int offset = warpSize / 2; offset > 0; offset /= 2) {
          total_trilinear += __shfl_down_sync(0xffffffff, total_trilinear, offset);
          total_eintcal += __shfl_down_sync(0xffffffff, total_eintcal, offset);
        }

        if (local_thread_id == 0) {
          const fp_type tors_free_energy        = num_rotamers * autodock_parameters::coeff_tors;
          s_chromosome_scores[chromosome_index] = total_trilinear + total_eintcal + tors_free_energy;
        }
      }

      // Generate the new population
      for (int chromosome_index = local_thread_id; chromosome_index < chromosome_number;
           chromosome_index += thread_per_block) {
        chromosome& next_chromosome = *(l_next_chromosomes + chromosome_index);

        // select the parent
        const int best_individual_1 =
            tournament_selection_cuda(l_state, tournament_length, chromosome_number, s_chromosome_scores);
        const int best_individual_2 =
            tournament_selection_cuda(l_state, tournament_length, chromosome_number, s_chromosome_scores);

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
      __syncwarp();
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
      const fp_type other_min_score = __shfl_down_sync(0xFFFFFFFF, min_score, offset);
      const int other_min_index     = __shfl_down_sync(0xFFFFFFFF, min_index, offset);
      if (other_min_score < min_score) {
        min_score = other_min_score;
        min_index = other_min_index;
      }
    }
    if (local_thread_id == 0) {
      ligand_scores[ligand_id] = min_score;
      memcpy((*(best_chromosomes + ligand_id)).data(),
             (*(l_chromosomes + min_index)).data(),
             sizeof(fp_type) * (6 + num_rotamers));
    }
  }
} // namespace mudock
