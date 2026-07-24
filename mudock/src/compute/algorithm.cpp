#include <algorithm>
#include <cctype>
#include <mudock/compute/algorithm.hpp>
#include <mudock/compute/pipeline_selector.hpp>
#include <string>

namespace mudock {
  search_algorithm parse_search_algorithm(const std::string_view token) {
    const auto element_it = std::find_if(std::begin(SEARCH_ALGORITHM_DICT),
                                         std::end(SEARCH_ALGORITHM_DICT),
                                         [&token](const auto& e) { return e.token == token; });
    if (element_it != std::end(SEARCH_ALGORITHM_DICT))
      return element_it->value;
    else
      throw std::runtime_error("Missing search algorithm");
  }
  scoring_function parse_scoring_function(const std::string_view token) {
    const auto element_it = std::find_if(std::begin(SCORING_FUNCTION_DICT),
                                         std::end(SCORING_FUNCTION_DICT),
                                         [&token](const auto& e) { return e.token == token; });
    if (element_it != std::end(SCORING_FUNCTION_DICT))
      return element_it->value;
    else
      throw std::runtime_error("Missing scoring function");
  }

  std::string_view to_string(const search_algorithm algorithm) {
    return SEARCH_ALGORITHM_DICT[static_cast<int>(algorithm)].token;
  }

  std::string_view to_string(const scoring_function scoring) {
    return SCORING_FUNCTION_DICT[static_cast<int>(scoring)].token;
  }
} // namespace mudock
