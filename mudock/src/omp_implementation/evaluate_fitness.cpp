#include <cstdio>
#include <mudock/grid.hpp>
#include <mudock/omp_implementation/calc_energy.hpp>
#include <mudock/omp_implementation/evaluate_fitness.hpp>
#include <mudock/omp_implementation/mutate.hpp>
#include <mudock/omp_implementation/omp_random.hpp>
#include <mudock/utils.hpp>
#include <omp.h>

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
  void evaluate_fitness(const int batch_ligands,
                        const int num_generations,
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
                        // const int* __restrict__ atom_tex_indexes,
                        // const cudaTextureObject_t electro_texture,
                        // const cudaTextureObject_t desolv_texture,
                        fp_type* __restrict__ ligand_scores,
                        chromosome* __restrict__ best_chromosomes) {
#pragma omp target teams num_teams(batch_ligands)
    for (auto ligand_id = 0; ligand_id < static_cast<int>(batch_ligands); ++ligand_id) {
      const int num_atoms    = ligand_num_atoms[ligand_id];
      const int num_nonbonds = ligand_num_nonbonds[ligand_id];
      const int num_rotamers = ligand_num_rotamers[ligand_id];
      int team_id            = omp_get_team_num();
      int num_teams          = omp_get_num_teams();
      int thread_id          = omp_get_thread_num();
      int num_threads        = omp_get_num_threads();

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
      // const auto* l_atom_tex_indexes      = atom_  tex_indexes + ligand_id * atom_stride;
      const int* l_ligand_nonbond_a1      = ligand_nonbond_a1 + ligand_id * nonbond_stride;
      const int* l_ligand_nonbond_a2      = ligand_nonbond_a2 + ligand_id * nonbond_stride;
      // XORWOWState l_state{};

#pragma omp parallel for
      for (int atom_index = 0; atom_index < num_atoms; ++atom_index) {
        // printf("%d | %d %d | %f %f %f\n",
        //        omp_get_thread_num(),
        //        ligand_id,
        //        atom_index,
        //        d_original_ligand_x[atom_index],
        //        d_original_ligand_y[atom_index],
        //        d_original_ligand_z[atom_index]);
      }
    }
  }
} // namespace mudock
