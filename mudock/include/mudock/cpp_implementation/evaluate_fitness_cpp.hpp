#pragma once

#include <cstring>
#include <memory>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/cpp_implementation/calc_energy.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/likwid_utils.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <random>
#include <stdint.h>

namespace mudock {
  static constexpr auto coordinate_step = fp_type{0.2};
  static constexpr auto angle_step      = fp_type{4};

  template<typename T>
  [[nodiscard]] inline const T random_gen_cpp(std::mt19937& generator,
                                              std::uniform_real_distribution<fp_type>& dist,
                                              const T& min,
                                              const T& max) {
    fp_type value;
    if constexpr (is_debug())
      value = fp_type{0.4};
    else {
      value = dist(generator);
    }
    return static_cast<T>(value * (max - min) + min);
  }

  inline int get_selection_distribution(std::mt19937& generator,
                                        std::uniform_real_distribution<fp_type>& dist,
                                        const int& population_number) {
    return random_gen_cpp<int>(generator, dist, 0, population_number - 1);
  };
  inline fp_type get_init_change_distribution(std::mt19937& generator,
                                              std::uniform_real_distribution<fp_type>& dist) {
    return random_gen_cpp<fp_type>(generator, dist, -45, 45);
  }
  inline fp_type get_mutation_change_distribution(std::mt19937& generator,
                                                  std::uniform_real_distribution<fp_type>& dist) {
    return random_gen_cpp<fp_type>(generator, dist, -10, 10);
  };
  inline fp_type get_mutation_coin_distribution(std::mt19937& generator,
                                                std::uniform_real_distribution<fp_type>& dist) {
    return random_gen_cpp<fp_type>(generator, dist, 0, 1);
  };
  inline int get_crossover_distribution(std::mt19937& generator,
                                        std::uniform_real_distribution<fp_type>& dist,
                                        const int& num_rotamers) {
    return random_gen_cpp<int>(generator, dist, 0, 6 + num_rotamers);
  };
  inline const chromosome& tournament_selection(std::mt19937& generator,
                                                std::uniform_real_distribution<fp_type>& dist,
                                                const int& tournament_length,
                                                const individual* __restrict__ population,
                                                const int& population_number) {
    const auto num_iterations = tournament_length;
    auto best_individual      = get_selection_distribution(generator, dist, population_number);
    for (int i = 0; i < num_iterations; ++i) {
      auto contendent = get_selection_distribution(generator, dist, population_number);
      if (population[contendent].score < population[best_individual].score) {
        best_individual = contendent;
      }
    }
    return population[best_individual].genes;
  }

  template<cpu_vectorization vect>
  void evaluate_fitness(const fp_type* __restrict__ ligand_x,
                        const fp_type* __restrict__ ligand_y,
                        const fp_type* __restrict__ ligand_z,
                        const fp_type* __restrict__ ligand_vol,
                        const fp_type* __restrict__ ligand_solpar,
                        const fp_type* __restrict__ ligand_charge,
                        const int* __restrict__ map_ligand_offsets,
                        const int num_atoms,
                        const int num_rotamers,
                        const int* __restrict__ frag_masks,
                        const int* __restrict__ frag_start_indexes,
                        const int* __restrict__ frag_stop_indexes,
                        const int num_nonbond,
                        const int* __restrict__ non_bond_list_a1,
                        const int* __restrict__ non_bond_list_a2,
                        const fp_type* __restrict__ cA_list,
                        const fp_type* __restrict__ cB_list,
                        const int* __restrict__ xB_list,
                        const fp_type* __restrict__ grid_maps,
                        const int num_generations,
                        const int population_size,
                        const int tournament_length,
                        const fp_type mutation_prob,
                        const fp_type* __restrict__ minimum,
                        const fp_type* __restrict__ maximum,
                        const fp_type* __restrict__ center,
                        const int map_index_x,
                        const int map_index_xy,
                        const int map_index_xyz,
                        individual* __restrict__ population_buffer1,
                        individual* __restrict__ population_buffer2,
                        const int seed) {
    std::uniform_real_distribution<fp_type> dist{fp_type{0.0}, fp_type{1.0}};
    std::mt19937 generator(seed);

    auto* population      = population_buffer1;
    auto* next_population = population_buffer2;
    std::array<fp_type, max_static_atoms()> altered_x{};
    std::array<fp_type, max_static_atoms()> altered_y{};
    std::array<fp_type, max_static_atoms()> altered_z{};

    // Randomly initialize the population
    for (int element_index = 0; element_index < population_size; ++element_index) {
      auto& element = population[element_index];
      for (int i{0}; i < 3; ++i) { // initialize the rigid translation
        element.genes[i] = get_init_change_distribution(generator, dist) * coordinate_step;
      }
      for (int i{3}; i < 6 + num_rotamers; ++i) { // initialize the rotations
        element.genes[i] = get_init_change_distribution(generator, dist) * angle_step;
      }
    }

    for (int generation = 0; generation < num_generations; ++generation) {
      LIKWID_MARKER_START("GA");
      // Evaluate the fitness of the population
      for (int element_index = 0; element_index < population_size; ++element_index) {
        auto& element = population[element_index];

        std::memcpy(altered_x.data(), ligand_x, num_atoms * sizeof(fp_type));
        std::memcpy(altered_y.data(), ligand_y, num_atoms * sizeof(fp_type));
        std::memcpy(altered_z.data(), ligand_z, num_atoms * sizeof(fp_type));

        // TODO check it it makes sense -> print the MOL2
        // apply the transformation encoded in the element genes to the original ligand
        apply<vect>(altered_x.data(),
                    altered_y.data(),
                    altered_z.data(),
                    element.genes,
                    num_atoms,
                    num_rotamers,
                    frag_masks,
                    frag_start_indexes,
                    frag_stop_indexes);

        // compute the energy of the system
        const auto energy = calc_energy<vect>(altered_x.data(),
                                              altered_y.data(),
                                              altered_z.data(),
                                              ligand_vol,
                                              ligand_solpar,
                                              ligand_charge,
                                              map_ligand_offsets,
                                              num_atoms,
                                              num_rotamers,
                                              num_nonbond,
                                              non_bond_list_a1,
                                              non_bond_list_a2,
                                              cA_list,
                                              cB_list,
                                              xB_list,
                                              minimum,
                                              maximum,
                                              center,
                                              map_index_x,
                                              map_index_xy,
                                              map_index_xyz,
                                              grid_maps);
        element.score     = energy; // dummy implementation to test the genetic
      }
      LIKWID_MARKER_STOP("GA");

      // Generate the new population
      for (int element_index = 0; element_index < population_size; ++element_index) {
        auto& next_individual = next_population[element_index];
        // select the parent
        auto best_individual_1 = get_selection_distribution(generator, dist, population_size);
        auto best_individual_2 = get_selection_distribution(generator, dist, population_size);
        for (int i = 0; i < tournament_length; ++i) {
          const auto contendent_1 = get_selection_distribution(generator, dist, population_size);
          const auto contendent_2 = get_selection_distribution(generator, dist, population_size);
          if (population[contendent_1].score < population[best_individual_1].score) {
            best_individual_1 = contendent_1;
          }
          if (population[contendent_2].score < population[best_individual_2].score) {
            best_individual_2 = contendent_2;
          }
        }
        const auto& parent1 = population[best_individual_1].genes;
        const auto& parent2 = population[best_individual_2].genes;

        // generate the offspring
        const auto split_index = get_crossover_distribution(generator, dist, num_rotamers);
        std::copy(std::begin(parent1), std::begin(parent1) + split_index, std::begin(next_individual.genes));
        std::copy(std::begin(parent2) + split_index,
                  std::end(parent2),
                  std::begin(next_individual.genes) + split_index);
        next_individual.score = fp_type{0};

        // mutate the offspring
        for (int i{0}; i < 3; ++i) {
          if (get_mutation_coin_distribution(generator, dist) < mutation_prob)
            next_individual.genes[i] += get_mutation_change_distribution(generator, dist) * coordinate_step;
        }
        for (int i{3}; i < 6 + num_rotamers; ++i) {
          if (get_mutation_coin_distribution(generator, dist) < mutation_prob)
            next_individual.genes[i] += get_mutation_change_distribution(generator, dist) * angle_step;
        }
      }

      // swap the new population with the old one
      const auto temp = population;
      population      = next_population;
      next_population = temp;
    }
  }
} // namespace mudock
