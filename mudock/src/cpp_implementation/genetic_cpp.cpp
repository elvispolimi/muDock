#include "mudock/type_alias.hpp"

#include <cstring>
#include <mudock/compute/buffer.hpp>
#include <mudock/compute/devices_memory.hpp>
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
  static constexpr auto coordinate_step = static_cast<fp_type>(0.2);
  static constexpr auto angle_step      = static_cast<fp_type>(4);

  thread_local device_memory<std::mt19937> rand_device;

  template<typename T>
  [[nodiscard]] inline const T random_gen_cpp(std::mt19937& generator,
                                              std::uniform_real_distribution<fp_type>& dist,
                                              const T& min,
                                              const T& max) {
    fp_type value;
    if constexpr (is_debug())
      value = static_cast<fp_type>(0.4);
    else {
      value = dist(generator);
    }
    return static_cast<T>(value * static_cast<fp_type>(max - min) + static_cast<fp_type>(min));
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

  void initialize_impl(const int batch_ligands,
                       const int population_number,
                       const int seed,
                       chromosome* population,
                       int* __restrict__ num_rotamers_b,
                       fp_type* __restrict__ scores_b) {
    std::uniform_real_distribution<fp_type> dist{fp_type{0.0}, fp_type{1.0}};
    rand_device.init(seed);
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      const int num_rotamers                = num_rotamers_b[ligand_index];
      chromosome* __restrict__ population_l = population + population_number * ligand_index;
      fp_type* __restrict__ scores          = scores_b + population_number * ligand_index;
      for (int element_index = 0; element_index < population_number; ++element_index) {
        chromosome& element   = population_l[element_index];
        scores[element_index] = big_bound<fp_type>();
        for (int i{0}; i < 3; ++i) { // initialize the rigid translation
          element[i] = get_init_change_distribution(rand_device(), dist) * coordinate_step;
        }
        for (int i{3}; i < 6 + num_rotamers; ++i) { // initialize the rotations
          element[i] = get_init_change_distribution(rand_device(), dist) * angle_step;
        }
      }
    }
  }
  void iterate_impl(const int batch_ligands,
                    const int population_number,
                    const int tournament_length,
                    const fp_type mutation_prob,
                    chromosome* population,
                    chromosome* next_population,
                    int* __restrict__ num_rotamers_b,
                    fp_type* __restrict__ scores_b) {
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
        auto best_individual_1 = get_selection_distribution(rand_device(), dist, population_number);
        auto best_individual_2 = get_selection_distribution(rand_device(), dist, population_number);
        for (int i = 0; i < tournament_length; ++i) {
          const auto contendent_1 = get_selection_distribution(rand_device(), dist, population_number);
          const auto contendent_2 = get_selection_distribution(rand_device(), dist, population_number);
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
        const auto split_index = get_crossover_distribution(rand_device(), dist, num_rotamers);
        std::copy(std::begin(parent1), std::begin(parent1) + split_index, std::begin(next_individual));
        std::copy(std::begin(parent2) + split_index,
                  std::end(parent2),
                  std::begin(next_individual) + split_index);

        // mutate the offspring
        for (int i{0}; i < 3; ++i) {
          if (get_mutation_coin_distribution(rand_device(), dist) < mutation_prob)
            next_individual[i] += get_mutation_change_distribution(rand_device(), dist) * coordinate_step;
        }
        for (int i{3}; i < 6 + num_rotamers; ++i) {
          if (get_mutation_coin_distribution(rand_device(), dist) < mutation_prob)
            next_individual[i] += get_mutation_change_distribution(rand_device(), dist) * angle_step;
        }
      }
    }
  }

  void finalize_impl(const int batch_ligands,
                     const int population_number,
                     chromosome* population,
                     int* __restrict__ num_rotamers_b,
                     fp_type* __restrict__ scores_b,
                     fp_type* __restrict__ best_scores_b,
                     chromosome* __restrict__ best_chromosomes_b) {
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
      best_scores_b[ligand_index] = min_score;
      memcpy((best_chromosomes_b + ligand_index)->data(),
             (population_l + min_index)->data(),
             sizeof(fp_type) * (6 + num_rotamers));
    }
  }

  template<>
  void genetic_kernel<queue_cpp>::finalize() {
    q->invoke_kernel<finalize_region_name>(finalize_impl,
                                           batch_ligands,
                                           population_number,
                                           population,
                                           num_rotamers_b,
                                           scores_b,
                                           best_scores_b,
                                           best_chromosomes_b);
  }
  template<>
  void genetic_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<iterate_region_name>(iterate_impl,
                                          batch_ligands,
                                          population_number,
                                          tournament_length,
                                          mutation_prob,
                                          population,
                                          next_population,
                                          num_rotamers_b,
                                          scores_b);
  }
  template<>
  void genetic_kernel<queue_cpp>::initialize() {
    q->invoke_kernel<initialize_region_name>(initialize_impl,
                                             batch_ligands,
                                             population_number,
                                             seed,
                                             population,
                                             num_rotamers_b,
                                             scores_b);
  }
} // namespace mudock
