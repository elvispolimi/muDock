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
  __device__ __constant__ fp_type inv_spacing_const;

  void setup_constant_memory(const point3D& minimum, const fp_type inv_grid_spacing) {
    fp_type l_map_min[3]{minimum.x, minimum.y, minimum.z};
    MUDOCK_CHECK(
        cudaMemcpyToSymbol(map_min_const, &l_map_min, 3 * sizeof(fp_type), 0, cudaMemcpyHostToDevice));
    MUDOCK_CHECK(
        cudaMemcpyToSymbol(inv_spacing_const, &inv_grid_spacing, sizeof(fp_type), 0, cudaMemcpyHostToDevice));
  }

  __device__ fp_type trilinear_interpolation_cuda(const fp_type coord[], const cudaTextureObject_t& tex) {
    // Interpolation CUDA
    fp_type p0u, p0v, p0w;

    const int u0      = (int) (coord[0]) + 1;
    const fp_type p1u = fp_type{1} - (p0u = coord[0] - static_cast<fp_type>(u0));

    const int v0      = (int) (coord[1]) + 1;
    const fp_type p1v = fp_type{1} - (p0v = coord[1] - static_cast<fp_type>(v0));

    const int w0      = (int) (coord[2]) + 1;
    const fp_type p1w = fp_type{1} - (p0w = coord[2] - static_cast<fp_type>(w0));

    const fp_type pu[2] = {p1u, p0u};
    const fp_type pv[2] = {p1v, p0v};
    const fp_type pw[2] = {p1w, p0w};
    fp_type value{0};
    for (int i = 0; i <= 1; i++)
      for (int t = 0; t <= 1; t++)
        for (int n = 0; n <= 1; n++) {
          const fp_type tmp = tex3D<fp_type>(tex, u0 + n, v0 + t, w0 + i);
          value += pu[n] * pv[t] * pw[i] * tmp;
        }
    return value;
  }

  template<typename T>
  __device__ const T random_gen_cuda(curandState* state, const T min, const T max) {
    T value;
    if constexpr (is_debug())
      // TODO value here for debug
      value = T{0.5};
    else {
      value = curand_uniform(state);
    }
    return static_cast<T>(value * (max - min) + min);
  }

  __device__ int get_selection_distribution(curandState* state, const int* population_number) {
    return random_gen_cuda<int>(state, 0, *population_number - 1);
  };

  __device__ int get_init_change_distribution(curandState* state) {
    return random_gen_cuda<int>(state, -45, 45);
  }
  __device__ int get_mutation_change_distribution(curandState* state) {
    return random_gen_cuda<int>(state, -10, 10);
  };
  __device__ fp_type get_mutation_coin_distribution(curandState* state) {
    return random_gen_cuda<fp_type>(state, 0, 1);
  };
  __device__ int get_crossover_distribution(curandState* state, const int* num_rotamers) {
    return random_gen_cuda<int>(state, 0, 6 + *num_rotamers);
  };

  __device__ chromosome* tournament_selection_cuda(curandState* state,
                                                   const int tournament_length,
                                                   chromosome* __restrict__ chromosomes,
                                                   const int chromosome_number,
                                                   fp_type* scores) {
    const int num_iterations = tournament_length;
    auto best_individual     = get_selection_distribution(state, &chromosome_number);
    for (std::size_t i = 0; i < num_iterations; ++i) {
      auto contendent = get_selection_distribution(state, &chromosome_number);
      if (scores[contendent] < scores[best_individual]) {
        best_individual = contendent;
      }
    }
    return chromosomes + best_individual;
  }

  //TODO template parameter based on number of atoms
  __global__ void evaluate_fitness(const int num_generations,
                                   const int tournament_length,
                                   const int mutation_prob,
                                   const int chromosome_number,
                                   const int chromosome_stride,
                                   const int atom_stride,
                                   const int rotamers_stride,
                                   const int nonbond_stride,
                                   fp_type* __restrict__ ligand_x,
                                   fp_type* __restrict__ ligand_y,
                                   fp_type* __restrict__ ligand_z,
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
                                   const cudaTextureObject_t* atom_textures,
                                   const int* atom_tex_indexes,
                                   const cudaTextureObject_t electro_texture,
                                   const cudaTextureObject_t desolv_texture,
                                   curandState* state,
                                   fp_type* ligand_scores) {
    const int ligand_id = blockIdx.x;
    const int threadId  = threadIdx.x + blockDim.x * blockIdx.x;

    const int num_atoms    = ligand_num_atoms[ligand_id];
    const int num_nonbonds = ligand_num_nonbonds[ligand_id];
    const int num_rotamers = ligand_num_rotamers[ligand_id];

    fp_type* l_ligand_x              = ligand_x + chromosome_number * ligand_id * atom_stride;
    fp_type* l_ligand_y              = ligand_y + chromosome_number * ligand_id * atom_stride;
    fp_type* l_ligand_z              = ligand_z + chromosome_number * ligand_id * atom_stride;
    const fp_type* l_ligand_vol      = ligand_vol + ligand_id * atom_stride;
    const fp_type* l_ligand_solpar   = ligand_solpar + ligand_id * atom_stride;
    const fp_type* l_ligand_charge   = ligand_charge + ligand_id * atom_stride;
    const int* l_ligand_num_hbond    = ligand_num_hbond + ligand_id * atom_stride;
    const fp_type* l_ligand_Rij_hb   = ligand_Rij_hb + ligand_id * atom_stride;
    const fp_type* l_ligand_Rii      = ligand_Rii + ligand_id * atom_stride;
    const fp_type* l_ligand_epsij_hb = ligand_epsij_hb + ligand_id * atom_stride;
    const fp_type* l_ligand_epsii    = ligand_epsii + ligand_id * atom_stride;
    chromosome* l_chromosomes        = chromosomes + ligand_id * chromosome_stride;
    // Point to the next population buffer
    chromosome* l_next_chromosomes      = chromosomes + ligand_id * chromosome_stride + chromosome_number;
    const auto* l_fragments             = ligand_fragments + ligand_id * atom_stride;
    const auto* l_frag_start_atom_index = frag_start_atom_index + ligand_id * rotamers_stride;
    const auto* l_frag_stop_atom_index  = frag_stop_atom_index + ligand_id * rotamers_stride;
    const auto* l_atom_tex_indexes      = atom_tex_indexes + ligand_id * atom_stride;
    const int* l_ligand_nonbond_a1      = ligand_nonbond_a1 + ligand_id * nonbond_stride;
    const int* l_ligand_nonbond_a2      = ligand_nonbond_a2 + ligand_id * nonbond_stride;
    curandState* l_state                = state + threadId;

    // Shared memory
    extern __shared__ fp_type shared_data[];
    fp_type* s_chromosome_scores = shared_data;

    // Generate initial population
    for (int chromosome_index = threadIdx.x; chromosome_index < chromosome_number;
         chromosome_index += blockDim.x) {
      chromosome* chromo = l_chromosomes + chromosome_index;
#pragma unroll
      for (int i{0}; i < 3; ++i) { // initialize the rigid translation
        (*chromo)[i] = static_cast<fp_type>(get_init_change_distribution(l_state)) * coordinate_step;
      }
#pragma unroll
      for (int i{3}; i < 3 + num_rotamers; ++i) { // initialize the rotations
        (*chromo)[i] = static_cast<fp_type>(get_init_change_distribution(l_state)) * angle_step;
      }
    }

    // TODO maybe template parameter?
    for (int generation = 0; generation < num_generations; ++generation) {
      for (int chromosome_index = 0; chromosome_index < chromosome_number; ++chromosome_index) {
        // Modify coordinates
        // TODO verify first the CPP version
        // apply_cuda(l_ligand_x,
        //            l_ligand_y,
        //            l_ligand_z,
        //            l_individual->genes.data(),
        //            l_fragments,
        //            l_frag_start_atom_index,
        //            l_frag_stop_atom_index,
        //            num_rotamers,
        //            atom_stride,
        //            num_atoms);

        __syncwarp();

        // Calculate energy
        fp_type elect_total_trilinear = 0;
        fp_type emap_total_trilinear  = 0;
        fp_type dmap_total_trilinear  = 0;
        for (int atom_index = threadIdx.x; atom_index < num_atoms; atom_index += blockDim.x) {
          fp_type coord_tex[3]{(l_ligand_x[atom_index] - map_min_const[0]) * inv_spacing_const,
                               (l_ligand_y[atom_index] - map_min_const[1]) * inv_spacing_const,
                               (l_ligand_z[atom_index] - map_min_const[2]) * inv_spacing_const};

          //  TODO check approximations with in hardware interpolation
          // elect_total_trilinear += tex3D<fp_type>(electro_texture, coord_tex[0], coord_tex[1], coord_tex[2]) *
          //                          l_ligand_charge[atom_index];
          // dmap_total_trilinear += tex3D<fp_type>(desolv_texture, coord_tex[0], coord_tex[1], coord_tex[2]);
          // emap_total_trilinear += tex3D<fp_type>(atom_textures[l_atom_tex_indexes[atom_index]],
          //                                        coord_tex[0],
          //                                        coord_tex[1],
          //                                        coord_tex[2]) *
          //                         fabsf(l_ligand_charge[atom_index]);

          elect_total_trilinear +=
              trilinear_interpolation_cuda(coord_tex, electro_texture) * l_ligand_charge[atom_index];
          dmap_total_trilinear +=
              trilinear_interpolation_cuda(coord_tex, desolv_texture) * fabsf(l_ligand_charge[atom_index]);
          emap_total_trilinear +=
              trilinear_interpolation_cuda(coord_tex, atom_textures[l_atom_tex_indexes[atom_index]]);
        }

        __syncwarp();

        fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
        if (num_rotamers > 0)
          calc_intra_energy(l_ligand_x,
                            l_ligand_y,
                            l_ligand_z,
                            l_ligand_vol,
                            l_ligand_solpar,
                            l_ligand_charge,
                            l_ligand_num_hbond,
                            l_ligand_Rij_hb,
                            l_ligand_Rii,
                            l_ligand_epsij_hb,
                            l_ligand_epsii,
                            num_atoms,
                            num_nonbonds,
                            l_ligand_nonbond_a1,
                            l_ligand_nonbond_a2,
                            &elect_total_eintcal,
                            &emap_total_eintcal,
                            &dmap_total_eintcal);

        __syncwarp();

        // Perform a tree reduction using __shfl_down_sync
        // TODO check performance
        for (int offset = 16; offset > 0; offset /= 2) {
          elect_total_trilinear += __shfl_down_sync(0xffffffff, elect_total_trilinear, offset);
          dmap_total_trilinear += __shfl_down_sync(0xffffffff, dmap_total_trilinear, offset);
          emap_total_trilinear += __shfl_down_sync(0xffffffff, emap_total_trilinear, offset);

          elect_total_eintcal += __shfl_down_sync(0xffffffff, elect_total_eintcal, offset);
          emap_total_eintcal += __shfl_down_sync(0xffffffff, emap_total_eintcal, offset);
          dmap_total_eintcal += __shfl_down_sync(0xffffffff, dmap_total_eintcal, offset);
        }
        if (threadIdx.x == 0) {
          const fp_type tors_free_energy        = num_rotamers * autodock_parameters::coeff_tors;
          s_chromosome_scores[chromosome_index] = elect_total_trilinear + dmap_total_trilinear +
                                                  emap_total_trilinear + elect_total_eintcal +
                                                  emap_total_eintcal + dmap_total_eintcal + tors_free_energy;
        }

        __syncwarp();
      }

      // Generate the new population
      for (int chromosome_index = threadIdx.x; chromosome_index < chromosome_number;
           chromosome_index += blockDim.x) {
        chromosome* next_chromosome = l_next_chromosomes + chromosome_index;
        // for (auto& next_individual: next_population) {
        // select the parent
        const chromosome* parent1 = tournament_selection_cuda(l_state,
                                                              tournament_length,
                                                              l_chromosomes,
                                                              chromosome_number,
                                                              s_chromosome_scores);
        const chromosome* parent2 = tournament_selection_cuda(l_state,
                                                              tournament_length,
                                                              l_chromosomes,
                                                              chromosome_number,
                                                              s_chromosome_scores);

        // generate the offspring
        const int split_index = get_crossover_distribution(l_state, &num_rotamers);
        std::memcpy(next_chromosome, parent1, split_index);
        std::memcpy(next_chromosome + split_index, parent2 + split_index, parent2->size() - split_index);

// mutate the offspring
#pragma unroll
        for (int i{0}; i < 3; ++i) {
          if (get_mutation_coin_distribution(l_state) < mutation_prob)
            (*next_chromosome)[i] +=
                static_cast<fp_type>(get_mutation_change_distribution(l_state)) * coordinate_step;
        }
#pragma unroll
        for (int i{3}; i < 3 + num_rotamers; ++i) {
          if (get_mutation_coin_distribution(l_state) < mutation_prob)
            (*next_chromosome)[i] +=
                static_cast<fp_type>(get_mutation_change_distribution(l_state)) * angle_step;
        }
      }

      // Swap the actual population with the next one
      chromosome* tmp_chromosomes = l_chromosomes;
      l_chromosomes               = l_next_chromosomes;
      l_next_chromosomes          = tmp_chromosomes;
    }
  }
} // namespace mudock
