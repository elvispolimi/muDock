#include <cstdio>
#include <mudock/grid.hpp>
#include <mudock/omp_implementation/calc_energy.hpp>
#include <mudock/omp_implementation/evaluate_fitness.hpp>
#include <mudock/omp_implementation/mutate.hpp>
#include <mudock/utils.hpp>

#define FLATTENED_3D(x, y, z, index) (index.size_xy() * z + y * index.size_x() + x)

namespace mudock {
  static constexpr fp_type coordinate_step{0.2};
  static constexpr fp_type angle_step{4};

  fp_type trilinear_interpolation_omp(const fp_type coord[], const fp_type* tex, const index3D& index) {
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
          const fp_type tmp = tex[FLATTENED_3D(u0 + n, v0 + t, w0 + i, index)];
          value += pu[n] * pv[t] * pw[i] * tmp;
        }
    return value;
  }

  // TODO random gen
  // template<typename T>
  // const T random_gen_cuda(curandState& state, const T min, const T max) {
  //   fp_type value;
  //   if constexpr (is_debug()) {
  //     // TODO value here for debug
  //     value = fp_type{0.4};
  //   } else {
  //     value = curand_uniform(&state);
  //   }
  //   return static_cast<T>((value * static_cast<fp_type>(max - min)) + min);
  // }

  // int get_selection_distribution(curandState& state, const int* population_number) {
  //   return random_gen_cuda<int>(state, 0, *population_number - 1);
  // };

  // fp_type get_init_change_distribution(curandState& state) {
  //   return random_gen_cuda<fp_type>(state, -45, 45);
  // }
  // fp_type get_mutation_change_distribution(curandState& state) {
  //   return random_gen_cuda<fp_type>(state, -10, 10);
  // };
  // fp_type get_mutation_coin_distribution(curandState& state) {
  //   return random_gen_cuda<fp_type>(state, 0, 1);
  // };
  // int get_crossover_distribution(curandState& state, const int* num_rotamers) {
  //   return random_gen_cuda<int>(state, 0, 6 + *num_rotamers);
  // };

  // TODO
  // int tournament_selection_cuda(curandState& state,
  //                               const int tournament_length,
  //                               const int chromosome_number,
  //                               const fp_type* __restrict__ scores) {
  //   const int num_iterations = tournament_length;
  //   int best_individual      = get_selection_distribution(state, &chromosome_number);
  //   for (int i = 0; i < num_iterations; ++i) {
  //     auto contended = get_selection_distribution(state, &chromosome_number);
  //     if (scores[contended] < scores[best_individual]) {
  //       best_individual = contended;
  //     }
  //   }
  //   return best_individual;
  // }

  // TODO check the syncwarp
  //TODO OPT: template parameter based on number of atoms, rotamers, chromosomes and population
  // Interesting the usage of the bucketizer
  void evaluate_fitness(const int num_generations,
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
                        // TODO
                        // const cudaTextureObject_t* __restrict__ atom_textures,
                        const int* __restrict__ atom_tex_indexes,
                        // const cudaTextureObject_t electro_texture,
                        // const cudaTextureObject_t desolv_texture,
                        // curandState* __restrict__ state,
                        fp_type* __restrict__ ligand_scores,
                        chromosome* __restrict__ best_chromosomes) {
    // TODO
  }
} // namespace mudock
