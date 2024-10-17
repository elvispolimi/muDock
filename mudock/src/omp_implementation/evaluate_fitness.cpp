#include <cstdio>
#include <math.h>
#include <mudock/grid.hpp>
#include <mudock/omp_implementation/calc_energy.hpp>
#include <mudock/omp_implementation/evaluate_fitness.hpp>
#include <mudock/omp_implementation/mutate.hpp>
#include <mudock/omp_implementation/omp_random.hpp>
#include <mudock/utils.hpp>
#include <omp.h>

#define FLATTENED_3D(x, y, z, index) (index.size_xy() * z + y * index.size_x() + x)

typedef struct {
  mudock::fp_type value;
  int index;
} min_index_pair;

#pragma omp declare reduction(min_index:MinIndexPair : omp_out =                     \
                                  (omp_in.value < omp_out.value) ? omp_in : omp_out) \
    initializer(omp_priv = {INFINITY, -1})

namespace mudock {
  static constexpr fp_type coordinate_step{0.2};
  static constexpr fp_type angle_step{4};

  fp_type trilinear_interpolation_omp(const fp_type coord[], const fp_type* tex, const index3D& index) {
    // Interpolation OpenMP
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

  template<typename T>
  const T random_gen_omp(XORWOWState& state, const T min, const T max) {
    fp_type value;
    if constexpr (is_debug()) {
      // TODO value here for debug
      value = fp_type{0.4};
    } else {
      value = state.next();
    }
    return static_cast<T>((value * static_cast<fp_type>(max - min)) + min);
  }

  int get_selection_distribution(XORWOWState& state, const int* population_number) {
    return random_gen_omp<int>(state, 0, *population_number - 1);
  };

  fp_type get_init_change_distribution(XORWOWState& state) { return random_gen_omp<fp_type>(state, -45, 45); }
  fp_type get_mutation_change_distribution(XORWOWState& state) {
    return random_gen_omp<fp_type>(state, -10, 10);
  };
  fp_type get_mutation_coin_distribution(XORWOWState& state) { return random_gen_omp<fp_type>(state, 0, 1); };
  int get_crossover_distribution(XORWOWState& state, const int* num_rotamers) {
    return random_gen_omp<int>(state, 0, 6 + *num_rotamers);
  };

  int tournament_selection_omp(XORWOWState& state,
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
                        fp_type* __restrict__ scratch_chromosome,
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
                        const point3D minimum,
                        const point3D maximum,
                        const point3D center,
                        const index3D index,
                        const fp_type* const __restrict__* const __restrict__ atom_textures,
                        const int* __restrict__ atom_tex_indexes,
                        const fp_type* __restrict__ electro_texture,
                        const fp_type* __restrict__ desolv_texture,
                        XORWOWState* __restrict__ state,
                        fp_type* __restrict__ ligand_scores,
                        chromosome* __restrict__ best_chromosomes) {
#pragma omp target teams distribute
    for (auto ligand_id = 0; ligand_id < static_cast<int>(batch_ligands); ++ligand_id) {
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
      // XORWOWState l_state{};

      // Initialize scratchpads and generate initial population
      fp_type* l_scratch_chromosome = scratch_chromosome + ligand_id * chromosome_number;
#pragma omp parallel for
      for (int chromosome_index = 0; chromosome_index < chromosome_number; ++chromosome_index) {
        l_scratch_chromosome[chromosome_index] =
            std::numeric_limits<fp_type>::infinity(); // Set initial score value
        chromosome& chromo = *(l_chromosomes + chromosome_index);
#pragma unroll
        for (int i{0}; i < 3; ++i) { // initialize the rigid translation
          chromo[i] = get_init_change_distribution(state[chromosome_index]) * coordinate_step;
        }
#pragma unroll
        for (int i{3}; i < 6 + num_rotamers; ++i) { // initialize the rotations
          chromo[i] = get_init_change_distribution(state[chromosome_index]) * angle_step;
        }
      }

      for (int generation = 0; generation < num_generations; ++generation) {
        for (int chromosome_index = 0; chromosome_index < chromosome_number; ++chromosome_index) {
#pragma omp parallel for
          for (int atom_index = 0; atom_index < num_atoms; ++atom_index) {
            l_scratch_ligand_x[atom_index] = l_original_ligand_x[atom_index];
            l_scratch_ligand_y[atom_index] = l_original_ligand_y[atom_index];
            l_scratch_ligand_z[atom_index] = l_original_ligand_z[atom_index];
          }

          // Modify coordinates
          // apply_cuda(l_scratch_ligand_x,
          //            l_scratch_ligand_y,
          //            l_scratch_ligand_z,
          //            *(l_chromosomes + chromosome_index),
          //            l_fragments,
          //            l_frag_start_atom_index,
          //            l_frag_stop_atom_index,
          //            num_rotamers,
          //            atom_stride,
          //            num_atoms);

          // Calculate energy
          fp_type elect_total_trilinear = 0;
          fp_type emap_total_trilinear  = 0;
          fp_type dmap_total_trilinear  = 0;
#pragma omp parallel for reduction(+ : elect_total_trilinear, emap_total_trilinear, dmap_total_trilinear)
          for (int atom_index = 0; atom_index < num_atoms; ++atom_index) {
            fp_type coord_tex[3]{l_scratch_ligand_x[atom_index],
                                 l_scratch_ligand_y[atom_index],
                                 l_scratch_ligand_z[atom_index]};

            if (coord_tex[0] < minimum.x || coord_tex[0] > maximum.x || coord_tex[1] < minimum.y ||
                coord_tex[1] > maximum.y || coord_tex[2] < minimum.z || coord_tex[2] > maximum.z) {
              // Is outside
              const fp_type distance_two = pow(fabs(coord_tex[0] - center.x), fp_type{2}) +
                                           pow(fabs(coord_tex[1] - center.y), fp_type{2}) +
                                           pow(fabs(coord_tex[2] - center.z), fp_type{2});

              const fp_type epenalty = distance_two * ENERGYPENALTY;
              elect_total_trilinear += epenalty;
              emap_total_trilinear += epenalty;
            } else {
              // Is inside
              // Center atom coordinates on the grid center
              coord_tex[0] = (l_scratch_ligand_x[atom_index] - minimum.x) * inv_spacing,
              coord_tex[1] = (l_scratch_ligand_y[atom_index] - minimum.y) * inv_spacing;
              coord_tex[2] = (l_scratch_ligand_z[atom_index] - minimum.z) * inv_spacing;
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

              elect_total_trilinear += trilinear_interpolation_omp(coord_tex, electro_texture, index) *
                                       l_ligand_charge[atom_index];
              dmap_total_trilinear += trilinear_interpolation_omp(coord_tex, desolv_texture, index) *
                                      fabsf(l_ligand_charge[atom_index]);
              emap_total_trilinear +=
                  trilinear_interpolation_omp(coord_tex,
                                              atom_textures[l_atom_tex_indexes[atom_index]],
                                              index);
            }
          }

          fp_type total_trilinear = elect_total_trilinear + dmap_total_trilinear + emap_total_trilinear;

          fp_type total_eintcal{0};
          // if (num_rotamers > 0)
          // total_eintcal += calc_intra_energy(l_scratch_ligand_x,
          //                                    l_scratch_ligand_y,
          //                                    l_scratch_ligand_z,
          //                                    l_ligand_vol,
          //                                    l_ligand_solpar,
          //                                    l_ligand_charge,
          //                                    l_ligand_num_hbond,
          //                                    l_ligand_Rij_hb,
          //                                    l_ligand_Rii,
          //                                    l_ligand_epsij_hb,
          //                                    l_ligand_epsii,
          //                                    num_nonbonds,
          //                                    l_ligand_nonbond_a1,
          //                                    l_ligand_nonbond_a2);

          const fp_type tors_free_energy         = num_rotamers * autodock_parameters::coeff_tors;
          l_scratch_chromosome[chromosome_index] = total_trilinear + total_eintcal + tors_free_energy;
        }

// Generate the new population
#pragma omp parallel for
        for (int chromosome_index = 0; chromosome_index < chromosome_number; ++chromosome_index) {
          chromosome& next_chromosome = *(l_next_chromosomes + chromosome_index);

          // select the parent
          const int best_individual_1 = tournament_selection_omp(state[chromosome_index],
                                                                  tournament_length,
                                                                  chromosome_number,
                                                                  l_scratch_chromosome);
          const int best_individual_2 = tournament_selection_omp(state[chromosome_index],
                                                                  tournament_length,
                                                                  chromosome_number,
                                                                  l_scratch_chromosome);

          // generate the offspring
          const int split_index = get_crossover_distribution(state[chromosome_index], &num_rotamers);
          memcpy(next_chromosome.data(),
                 &(l_chromosomes[best_individual_1][0]),
                 split_index * sizeof(fp_type));
          const int parent2_copy_size = 6 + num_rotamers - split_index;
          if (parent2_copy_size > 0)
            memcpy(next_chromosome.data() + split_index,
                   &(l_chromosomes[best_individual_2][split_index]),
                   parent2_copy_size * sizeof(fp_type));

// mutate the offspring
#pragma unroll
          for (int i{0}; i < 3; ++i) {
            if (get_mutation_coin_distribution(state[chromosome_index]) < mutation_prob)
              next_chromosome[i] +=
                  get_mutation_change_distribution(state[chromosome_index]) * coordinate_step;
          }
#pragma unroll
          for (int i{3}; i < 6 + num_rotamers; ++i) {
            if (get_mutation_coin_distribution(state[chromosome_index]) < mutation_prob) {
              next_chromosome[i] += get_mutation_change_distribution(state[chromosome_index]) * angle_step;
            }
          }
        }

        // Swap the actual population with the next one
        chromosome* const tmp_chromosomes = l_chromosomes;
        l_chromosomes                     = l_next_chromosomes;
        l_next_chromosomes                = tmp_chromosomes;
      }

      // Output
      // Parallel for with maximum reduction
      min_index_pair result = {std::numeric_limits<fp_type>::infinity(), -1};
#pragma omp parallel for reduction(min_index : result)
      for (int chromosome_index = 0; chromosome_index < chromosome_number; ++chromosome_index) {
        if (l_scratch_chromosome[chromosome_index] < result.value) {
          result.value = l_scratch_chromosome[chromosome_index];
          result.index = chromosome_index;
        }
      }

      ligand_scores[ligand_id] = result.value;
      memcpy((*(best_chromosomes + ligand_id)).data(),
             (*(l_chromosomes + result.index)).data(),
             sizeof(fp_type) * (6 + num_rotamers));
    }
  }
} // namespace mudock
