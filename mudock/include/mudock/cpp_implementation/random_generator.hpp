#pragma once

#include <array>
#include <cstddef>
#include <random>

namespace mudock {
  static constexpr std::size_t random_number_size{1024};

  template<typename T>
  struct random_generator {
    random_generator() = default;
    random_generator(const std::size_t num_atoms,
                     const std::size_t num_rotamers,
                     const std::size_t population_size) {
      std::mt19937 generator{num_atoms};

      auto _init_change_distribution = std::uniform_int_distribution(-45, 45);
      for (std::size_t index = 0; index < random_number_size; ++index)
        init_change_distribution[index] = _init_change_distribution(generator);

      auto _selection_distribution = std::uniform_int_distribution<T>(std::size_t{0}, population_size - 1);
      for (std::size_t index = 0; index < random_number_size; ++index)
        selection_distribution[index] = _selection_distribution(generator);

      auto _mutation_change_distribution = std::uniform_int_distribution(-10, 10);
      for (std::size_t index = 0; index < random_number_size; ++index)
        mutation_change_distribution[index] = _mutation_change_distribution(generator);

      auto _mutation_coin_distribution = std::uniform_real_distribution{fp_type{0}, fp_type{1.0}};
      for (std::size_t index = 0; index < random_number_size; ++index)
        mutation_coin_distribution[index] = _mutation_coin_distribution(generator);

      auto _crossover_distribution =
          std::uniform_int_distribution<std::size_t>(std::size_t{0}, std::size_t{6} + num_rotamers);
      for (std::size_t index = 0; index < random_number_size; ++index)
        crossover_distribution[index] = _crossover_distribution(generator);
    };

    [[nodicard]] inline const T& get_init_change_distribution() const {
      index_init_change_distribution = index_init_change_distribution++ % random_number_size;
      return init_change_distribution[index_init_change_distribution];
    };

    [[nodicard]] inline const T& get_selection_distribution() const {
      index_selection_distribution = index_selection_distribution++ % random_number_size;
      return selection_distribution[index_selection_distribution];
    };
    [[nodicard]] inline const T& get_mutation_change_distribution() const {
      index_mutation_change_distribution = index_mutation_change_distribution++ % random_number_size;
      return mutation_change_distribution[index_mutation_change_distribution];
    };
    [[nodicard]] inline const T& get_mutation_coin_distribution() const {
      index_mutation_coin_distribution = index_mutation_coin_distribution++ % random_number_size;
      return mutation_coin_distribution[index_mutation_coin_distribution];
    };
    [[nodicard]] inline const T& get_crossover_distribution() const {
      index_crossover_distribution = index_crossover_distribution++ % random_number_size;
      return crossover_distribution[index_crossover_distribution];
    };

  private:
    std::array<T, random_number_size> init_change_distribution;
    std::array<T, random_number_size> selection_distribution;
    std::array<T, random_number_size> mutation_change_distribution;
    std::array<T, random_number_size> mutation_coin_distribution;
    std::array<T, random_number_size> crossover_distribution;
    int index_init_change_distribution{0}, index_selection_distribution{0},
        index_mutation_change_distribution{0}, index_mutation_coin_distribution{0},
        index_crossover_distribution{0};
  };
} // namespace mudock