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

  // template <typename fp_type>
  // std::vector<int>
  // get_indices_of_n_best(const fp_type* scores, int population_size, int n_best){
  //
  //   // Indices [0, 1, ..., population_size-1]
  //   std::vector<int> indices(population_size);
  //   std::iota(indices.begin(), indices.end(), 0);
  // 
  //   // Put the n_best lowest scores at the beginning.
  //   std::nth_element(
  //     indices.begin(),
  //     indices.begin() + n_best,
  //     indices.end(),
  //     [scores](int a, int b) {
  //         return scores[a] < scores[b];
  //     }
  //   );
  // 
  //   indices.resize(n_best);
  // 
  //   // Keep only the best 10%.
  //   indices.resize(n_best);
  // 
  //   return {std::move(indices)};
  // }
  // 
  // fp_type calculate_rmsd(chromosome first_individual, 
  //                        chromosome second_individual,
  //                        const int num_atoms) {
  //   fp_type rmsd = 0.0;
  //   for (int i = 0; i < num_atoms; ++i) {
  //     fp_type d2 = pow((first_individual[0]-second_individual[0]), 2) +
  //                 pow((first_individual[1]-second_individual[1]), 2) +
  //                 pow((first_individual[2]-second_individual[2]), 2)
  //   }
  // }

  // TODO L 
  // 1. it's not considering the crystal
  // 2. best_score_improved has a window of 1, maybe too rigid
  // 3. not sure about population_is_diverse criterion, seems too stupid for convergence
  bool has_converged(const fp_type best_so_far, 
                     const fp_type this_gen_best,
                     const fp_type best_score_diff_thld,
                     fp_type* __restrict__ scores,
                     const int population_size,
                     const fp_type score_variance_thld) {
    // Criterion 1: best score so far has not improved
    bool best_score_not_improved = (best_so_far - this_gen_best < best_score_diff_thld) ? true : false;
    if (best_score_not_improved) printf("Convergence due to best score not improving\n");

    // Criterion 2: population score is diverse
    fp_type mean = 0.0;
    for (int i = 0; i < population_size; ++i)
      mean += scores[i];
    mean /= static_cast<fp_type>(population_size);

    fp_type var = 0.0;
    for (int i = 0; i < population_size; ++i) {
      const fp_type d = scores[i] - mean;
      var += d * d;
    }
    var /= static_cast<fp_type>(population_size);
    
    bool population_score_is_similar = (var < score_variance_thld) ? true : false;
    if (population_score_is_similar) printf("Convergence due to population score not diverse\n");

    // WARNING: didn't implement this criterion because i don't have the atoms' coordinates here (x/y/z_scratch)
    // Criterion 3: RMSD of best 10%
    // const int n_best = std::max(1, population_size / 10);
    // std::vector<int> indices = get_indices_of_n_best(scores, population_size, n_best);
    // int similars = 0;
    // for (int i : indices) {
    //   const fp_type rmsd = calculate_rmsd(population[best_index], population[i], num_atoms);
    //   if (rmsd < rmsd_thld) {
    //     similars++;
    //   }
    // }
    // const fp_type fraction =  static_cast<fp_type>(similar) / n_best;
    // bool rmsd_is_low = (fraction >= 0.8) ? true : false;
    // if (rmsd_is_low) printf("Convergence due to low rmsd\n");

    // If something is triggered -> converge
    return best_score_not_improved || population_score_is_similar;
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
                    const int elite_size,
                    const int tournament_length,
                    const fp_type mutation_prob,
                    chromosome* population,
                    chromosome* next_population,
                    int* __restrict__ num_rotamers_b,
                    fp_type* __restrict__ scores_b,
                    fp_type* __restrict__ best_so_far_b,
                    const int generation,
                    const fp_type score_variance_thld,
                    const fp_type best_score_diff_thld,
                    const bool autostop,
                    const fp_type crystal_score,
                    const fp_type crystal_tolerance,
                    int* __restrict__ converged_ligands_b) {
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      const int converged_ligand = converged_ligands_b[ligand_index];
      if (converged_ligand) {
        continue;
      }
      std::uniform_real_distribution<fp_type> dist{fp_type{0.0}, fp_type{1.0}};
      const int num_rotamers                     = num_rotamers_b[ligand_index];
      chromosome* __restrict__ population_l      = population + population_number * ligand_index;
      chromosome* __restrict__ next_population_l = next_population + population_number * ligand_index;
      fp_type* __restrict__ scores               = scores_b + population_number * ligand_index;
      const fp_type best_so_far                  = best_so_far_b[ligand_index];

      // print best score
      fp_type best = scores[0];
      for(int i = 0; i < population_number; ++i){
        if (scores[i] < best){
          best = scores[i];
        }
      }
      printf("Gen %d --- Best of this gen: %f --- Best so far: %f\n", generation, double(best), double(best_so_far));
      // end print best score 


      // TODO L move autostop logic inside genetic.hpp maybe as separated stage, maybe doing it every n generations instead of doing it every generation
      // Check convergence 
      if (autostop && has_converged(best_so_far, best, best_score_diff_thld, scores, population_number, score_variance_thld)) {
        converged_ligands_b[ligand_index] = generation;
      }

      // Update best so far
      if (best < best_so_far) {
        best_so_far_b[ligand_index] = best;
      }

      // Elitism: preserve the best elite_size individuals
      // elite_indices[k] contains the index in population_l of the k-th best individual
      std::vector<int> elite_indices(elite_size, -1);

      for (int i = 0; i < population_number; ++i) {
        for (int e = 0; e < elite_size; ++e) {
          if (elite_indices[e] == -1 || scores[i] < scores[elite_indices[e]]) {

            // shift worse elites to the right
            for (int shift = elite_size - 1; shift > e; --shift) {
              elite_indices[shift] = elite_indices[shift - 1];
            }

            elite_indices[e] = i;
            break;
          }
        }
      }

      for (int e = 0; e < elite_size; ++e) {
        if (elite_indices[e] == -1)
          break;

        std::copy(std::begin(population_l[elite_indices[e]]),
                  std::end(population_l[elite_indices[e]]),
                  std::begin(next_population_l[e]));
      }

      // Generate the new population
      for (int element_index = elite_size; element_index < population_number; ++element_index) {
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
    q->invoke_kernel<this->iterate_region_name>(iterate_impl,
                                                batch_ligands,
                                                population_number,
                                                elite_size,
                                                tournament_length,
                                                mutation_prob,
                                                population,
                                                next_population,
                                                num_rotamers_b,
                                                scores_b,
                                                best_so_far_b,
                                                current_generation,
                                                score_variance_thld,
                                                best_score_diff_thld,
                                                autostop,
                                                crystal_score,
                                                crystal_tolerance,
                                                converged_ligands_b);
    ++current_generation;
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
