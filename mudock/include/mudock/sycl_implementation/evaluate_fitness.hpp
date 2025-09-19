#pragma once

#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/grid.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/sycl_implementation/calc_energy.hpp>
#include <mudock/sycl_implementation/mutate.hpp>
#include <mudock/sycl_implementation/sycl_random.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  // TODO move to a coeff file
  static constexpr fp_type coordinate_step{0.2};
  static constexpr fp_type angle_step{4};

  // TODO check the real randomness
  template<typename T>
  const T random_gen_sycl(XORWOWState& state, const T min, const T max) {
    fp_type value;
    if constexpr (is_debug())
      // TODO value here for debug
      value = fp_type{0.4};
    else {
      value = state.next();
    }
    return static_cast<T>((value * static_cast<fp_type>(max - min)) + min);
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

  template<int MAX_ATOMS>
  struct evaluate_fitness_kernel_tag {}; // just a name, no members

  template<int MAX_ATOMS>
  void evaluate_fitness(const int num_generations,
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
                        const point3D minimum,
                        const point3D maximum,
                        const point3D center,
                        const fp_type* __restrict__ grid_maps,
                        const int* __restrict__ map_ligand_offsets,
                        XORWOWState* __restrict__ state,
                        sycl::local_accessor<fp_type> s_chromosome_scores_acc,
                        // fp_type* __restrict__ s_chromosome_scores,
                        fp_type* __restrict__ ligand_scores,
                        chromosome* __restrict__ best_chromosomes,
                        sycl::nd_item<1> it) {
    const int workgroup_size       = it.get_local_range(0);
    const int workgroup_id         = it.get_group(0);
    const int ligand_id            = workgroup_id;
    const int workitem_id_in_group = it.get_local_id(0);
    const int workitem_id          = workitem_id_in_group + workgroup_size * workgroup_id;
    const auto& sub_group          = it.get_sub_group();

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

    fp_type* s_chromosome_scores = s_chromosome_scores_acc.get_multi_ptr<sycl::access::decorated::no>().get();

    XORWOWState& l_state = (state[workitem_id]);

    // Initialize score
    // TODO check with different workitem sizes
    for (int chromosome_index = workitem_id_in_group;
         chromosome_index < sycl::max(chromosome_number, workgroup_size);
         chromosome_index += workgroup_size)
      // FIXME undefined behaviour due to floating point flags
      s_chromosome_scores[chromosome_index] =
          std::numeric_limits<fp_type>::infinity(); // Set initial score value

    // Generate initial population
    for (int chromosome_index = workitem_id_in_group; chromosome_index < chromosome_number;
         chromosome_index += workgroup_size) {
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
    sycl::group_barrier(sub_group, sycl::memory_scope::sub_group);

    for (int generation = 0; generation < num_generations; ++generation) {
      for (int chromosome_index = 0; chromosome_index < chromosome_number; ++chromosome_index) {
        // Copy original coordinates
        // TODO shared memory?
#pragma unroll
        for (int atom_index = workitem_id_in_group; atom_index < MAX_ATOMS; atom_index += workgroup_size) {
          if (atom_index < num_atoms) {
            l_scratch_ligand_x[atom_index] = l_original_ligand_x[atom_index];
            l_scratch_ligand_y[atom_index] = l_original_ligand_y[atom_index];
            l_scratch_ligand_z[atom_index] = l_original_ligand_z[atom_index];
          }
        }
        // Modify coordinates
        apply_sycl<MAX_ATOMS>(l_scratch_ligand_x,
                              l_scratch_ligand_y,
                              l_scratch_ligand_z,
                              *(l_chromosomes + chromosome_index),
                              l_fragments,
                              l_frag_start_atom_index,
                              l_frag_stop_atom_index,
                              num_rotamers,
                              atom_stride,
                              num_atoms,
                              it);

        auto total_eintcal = calc_intra_energy<MAX_ATOMS>(l_scratch_ligand_x,
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
                                                          minimum,
                                                          maximum,
                                                          center,
                                                          map_index_x,
                                                          map_index_xy,
                                                          map_index_xyz,
                                                          grid_maps,
                                                          l_map_ligand_offsets,
                                                          it);

        total_eintcal = sycl::reduce_over_group(sub_group, total_eintcal, std::plus<fp_type>());
        if (workitem_id_in_group == 0) {
          // const fp_type tors_free_energy        = num_rotamers * autodock_parameters::coeff_tors;
          s_chromosome_scores[chromosome_index] = total_eintcal;
        }
      }

      // Generate the new population
      for (int chromosome_index = workitem_id_in_group; chromosome_index < chromosome_number;
           chromosome_index += workgroup_size) {
        chromosome& next_chromosome = *(l_next_chromosomes + chromosome_index);

        // select the parent
        const int best_individual_1 =
            tournament_selection_sycl(l_state, tournament_length, chromosome_number, s_chromosome_scores);
        const int best_individual_2 =
            tournament_selection_sycl(l_state, tournament_length, chromosome_number, s_chromosome_scores);

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
      sycl::group_barrier(sub_group, sycl::memory_scope::sub_group);
    }

    // Output

    // Compute the maximum value within the warp
    int min_index     = workitem_id_in_group;
    fp_type min_score = s_chromosome_scores[min_index];

    for (int chromosome_index = workitem_id_in_group + workgroup_size; chromosome_index < chromosome_number;
         chromosome_index += workgroup_size) {
      if (min_score > s_chromosome_scores[chromosome_index]) {
        min_index = chromosome_index;
        min_score = s_chromosome_scores[chromosome_index];
      }
    }

    const fp_type best_score = sycl::reduce_over_group(sub_group, min_score, sycl::minimum());

    // TODO checks that only one can do it
    if (min_score == best_score) {
      ligand_scores[ligand_id] = min_score + num_rotamers * autodock_parameters::coeff_tors;
      memcpy((*(best_chromosomes + ligand_id)).data(),
             (*(l_chromosomes + min_index)).data(),
             sizeof(fp_type) * (6 + num_rotamers));
    }
  };
} // namespace mudock
