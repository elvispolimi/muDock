#pragma once

#include <array>
#include <string_view>
namespace mudock {

  enum class search_algorithm : int { NONE = 0, GENETIC };
  enum class scoring_function : int { ADT = 0 };

  struct search_algorithm_description {
    search_algorithm value;
    std::string_view token;
  };

  struct scoring_function_description {
    scoring_function value;
    std::string_view token;
  };

  inline constexpr std::array<search_algorithm_description, 2> SEARCH_ALGORITHM_DICT = {
      {{search_algorithm::NONE, "none"}, {search_algorithm::GENETIC, "genetic"}}};

  inline constexpr std::array<scoring_function_description, 1> SCORING_FUNCTION_DICT = {
      {{scoring_function::ADT, "adt"}}};

  [[nodiscard]] search_algorithm parse_search_algorithm(std::string_view token);
  [[nodiscard]] scoring_function parse_scoring_function(std::string_view token);
  [[nodiscard]] std::string_view to_string(search_algorithm algorithm);
  [[nodiscard]] std::string_view to_string(scoring_function scoring);

  [[nodiscard]] constexpr int num_search_algorithms() {
    return static_cast<int>(SEARCH_ALGORITHM_DICT.size());
  }

  [[nodiscard]] constexpr int num_scoring_functions() {
    return static_cast<int>(SCORING_FUNCTION_DICT.size());
  }
} // namespace mudock
