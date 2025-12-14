#include <cstring>
#include <mudock/compute/buffer.hpp>
#include <mudock/compute/genetic.hpp>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cpp_implementation/genetic_cpp.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/likwid_utils.hpp>
#include <mudock/utils.hpp>
#include <random>

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

  template<>
  void genetic_kernel<queue_cpp>::initialize() {
    std::uniform_real_distribution<fp_type> dist{fp_type{0.0}, fp_type{1.0}};
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      const int num_rotamers                = num_rotamers_b[ligand_index];
      chromosome* __restrict__ population_l = population + population_number * ligand_index;
      fp_type* __restrict__ scores          = scores_b + population_number * ligand_index;
      for (int element_index = 0; element_index < population_number; ++element_index) {
        chromosome& element   = population_l[element_index];
        scores[element_index] = big_bound<fp_type>();
        for (int i{0}; i < 3; ++i) { // initialize the rigid translation
          element[i] = get_init_change_distribution(rand, dist) * coordinate_step;
        }
        for (int i{3}; i < 6 + num_rotamers; ++i) { // initialize the rotations
          element[i] = get_init_change_distribution(rand, dist) * angle_step;
        }
      }
    }
  }
  template<>
  void genetic_kernel<queue_cpp>::operator()() {
    // for (size_t generation = 0; generation < num_generations; ++generation) {
    //   LIKWID_MARKER_START("GA");
    //   // Evaluate the fitness of the population
    //   geom_trans();
    //   (*score_stage)();
    //
    //   LIKWID_MARKER_STOP("GA");

    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      std::uniform_real_distribution<fp_type> dist{fp_type{0.0}, fp_type{1.0}};
      const int num_rotamers                     = num_rotamers_b[ligand_index];
      chromosome* __restrict__ population_l      = population + population_number * ligand_index;
      chromosome* __restrict__ next_population_l = next_population + population_number * ligand_index;
      fp_type* __restrict__ scores               = scores_b + population_number * ligand_index;
      // Generate the new population
      for (int element_index = 0; element_index < population_number; ++element_index) {
        auto& next_individual = next_population_l[element_index];
        // select the parent
        auto best_individual_1 = get_selection_distribution(rand, dist, population_number);
        auto best_individual_2 = get_selection_distribution(rand, dist, population_number);
        for (int i = 0; i < tournament_length; ++i) {
          const auto contendent_1 = get_selection_distribution(rand, dist, population_number);
          const auto contendent_2 = get_selection_distribution(rand, dist, population_number);
          if (scores[contendent_1] < scores[best_individual_1]) {
            best_individual_1 = contendent_1;
          }
          if (scores[contendent_2] < scores[best_individual_2]) {
            best_individual_2 = contendent_2;
          }
        }
        const auto& parent1 = population_l[best_individual_1];
        const auto& parent2 = population_l[best_individual_2];

        // generate the offspring
        const auto split_index = get_crossover_distribution(rand, dist, num_rotamers);
        std::copy(std::begin(parent1), std::begin(parent1) + split_index, std::begin(next_individual));
        std::copy(std::begin(parent2) + split_index,
                  std::end(parent2),
                  std::begin(next_individual) + split_index);

        // mutate the offspring
        for (int i{0}; i < 3; ++i) {
          if (get_mutation_coin_distribution(rand, dist) < mutation_prob)
            next_individual[i] += get_mutation_change_distribution(rand, dist) * coordinate_step;
        }
        for (int i{3}; i < 6 + num_rotamers; ++i) {
          if (get_mutation_coin_distribution(rand, dist) < mutation_prob)
            next_individual[i] += get_mutation_change_distribution(rand, dist) * angle_step;
        }
      }
    }
  }
  template<>
  void genetic_kernel<queue_cpp>::finalize() {
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      const int num_rotamers                = num_rotamers_b[ligand_index];
      chromosome* __restrict__ population_l = population + population_number * ligand_index;
      fp_type* __restrict__ scores          = scores_b + population_number * ligand_index;
      int min_index                         = 0;
      fp_type min_score                     = scores[0];
      for (int chromosome_index = 0; chromosome_index < population_number; chromosome_index++) {
        if (min_score > scores[chromosome_index]) {
          min_index = chromosome_index;
          min_score = scores[chromosome_index];
        }
      }
      const fp_type tors_free_energy = num_rotamers * autodock_parameters::coeff_tors;
      best_scores_b[ligand_index]    = min_score + tors_free_energy;
      memcpy((best_chromosomes_b + ligand_index)->data(),
             (population_l + min_index)->data(),
             sizeof(fp_type) * (6 + num_rotamers));
    }
  }
} // namespace mudock
