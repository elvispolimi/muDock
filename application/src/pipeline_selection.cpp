#include "pipeline_selection.hpp"

#include <algorithm>
#include <cctype>
#include <string>

std::optional<search_algorithm> parse_search_algorithm(const std::string_view token) {
  std::string normalized{token};
  std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](const unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  const auto element_it =
      std::find_if(std::begin(SEARCH_ALGORITHM_TOKENS),
                   std::end(SEARCH_ALGORITHM_TOKENS),
                   [&normalized](const auto& element) { return element.token == normalized; });
  if (element_it != std::end(SEARCH_ALGORITHM_TOKENS)) {
    return element_it->value;
  }
  return std::nullopt;
}

std::optional<scoring_function> parse_scoring_function(const std::string_view token) {
  std::string normalized{token};
  std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](const unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  const auto element_it =
      std::find_if(std::begin(SCORING_FUNCTION_TOKENS),
                   std::end(SCORING_FUNCTION_TOKENS),
                   [&normalized](const auto& element) { return element.token == normalized; });
  if (element_it != std::end(SCORING_FUNCTION_TOKENS)) {
    return element_it->value;
  }
  return std::nullopt;
}

std::string_view to_string(const search_algorithm algorithm) {
  const auto element_it = std::find_if(std::begin(SEARCH_ALGORITHM_TOKENS),
                                       std::end(SEARCH_ALGORITHM_TOKENS),
                                       [&algorithm](const auto& element) { return element.value == algorithm; });
  if (element_it != std::end(SEARCH_ALGORITHM_TOKENS)) {
    return element_it->token;
  }
  return "unknown";
}

std::string_view to_string(const scoring_function scoring) {
  const auto element_it = std::find_if(std::begin(SCORING_FUNCTION_TOKENS),
                                       std::end(SCORING_FUNCTION_TOKENS),
                                       [&scoring](const auto& element) { return element.value == scoring; });
  if (element_it != std::end(SCORING_FUNCTION_TOKENS)) {
    return element_it->token;
  }
  return "unknown";
}
