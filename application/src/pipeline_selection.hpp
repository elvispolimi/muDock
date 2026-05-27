#pragma once

#include <mudock/compute/pipeline.hpp>
#include <mudock/utils.hpp>
#include <array>
#include <optional>
#include <string_view>

enum class search_algorithm : int { NONE = 0, GENETIC, COUNT };
enum class scoring_function : int { ADT = 0, COUNT };

struct search_algorithm_description {
  search_algorithm value;
  std::string_view token;
};

struct scoring_function_description {
  scoring_function value;
  std::string_view token;
};

static constexpr std::array<search_algorithm_description, 2> SEARCH_ALGORITHM_TOKENS = {
    {{search_algorithm::NONE, "none"}, {search_algorithm::GENETIC, "genetic"}}};

static constexpr std::array<scoring_function_description, 1> SCORING_FUNCTION_TOKENS = {
    {{scoring_function::ADT, "adt"}}};

[[nodiscard]] std::optional<search_algorithm> parse_search_algorithm(std::string_view token);
[[nodiscard]] std::optional<scoring_function> parse_scoring_function(std::string_view token);
[[nodiscard]] std::string_view to_string(search_algorithm algorithm);
[[nodiscard]] std::string_view to_string(scoring_function scoring);

[[nodiscard]] constexpr int num_search_algorithms() {
  return static_cast<int>(search_algorithm::COUNT);
}

[[nodiscard]] constexpr int num_scoring_functions() {
  return static_cast<int>(scoring_function::COUNT);
}

[[nodiscard]] constexpr int compose_pipeline_index(const search_algorithm search,
                                                   const scoring_function scoring) {
  return static_cast<int>(search) * num_scoring_functions() + static_cast<int>(scoring);
}

template<search_algorithm search, scoring_function scoring>
struct pipeline_selector;

template<>
struct pipeline_selector<search_algorithm::NONE, scoring_function::ADT> {
  using type = mudock::adt_score_pipeline;
};

template<>
struct pipeline_selector<search_algorithm::GENETIC, scoring_function::ADT> {
  using type = mudock::genetic_adt_pipeline;
};

template<search_algorithm search, scoring_function scoring>
using pipeline_selector_t = typename pipeline_selector<search, scoring>::type;

template<typename callback_t>
void dispatch_selected_pipeline(const search_algorithm search,
                                const scoring_function scoring,
                                callback_t&& callback) {
  constexpr_switch<0, num_search_algorithms(), 1>(
      [&](const auto search_index) {
        constexpr auto selected_search = static_cast<search_algorithm>(decltype(search_index)::value);
        constexpr_switch<0, num_scoring_functions(), 1>(
            [&](const auto scoring_index) {
              constexpr auto selected_scoring =
                  static_cast<scoring_function>(decltype(scoring_index)::value);
              using pipeline_t = pipeline_selector_t<selected_search, selected_scoring>;
              callback.template operator()<pipeline_t>(selected_search, selected_scoring);
            },
            static_cast<int>(scoring));
      },
      static_cast<int>(search));
}
