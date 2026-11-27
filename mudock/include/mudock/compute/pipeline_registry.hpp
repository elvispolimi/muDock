#pragma once

#include <mudock/compute/algorithm.hpp>
#include <mudock/compute/pipeline.hpp>
#include <tuple>

namespace mudock {
  template<typename... pipeline_types>
  struct pipeline_registry {};

  using scoring_pipeline_registry = pipeline_registry<scoring_pipeline<scoring_function::ADT>>;
  using genetic_pipeline_registry = pipeline_registry<genetic_scoring_pipeline<scoring_function::ADT>>;

  template<typename backend, scoring_function... scorers>
  auto make_genetic_backend_references(pipeline_registry<genetic_scoring_pipeline<scorers>...>) {
    return std::tuple{&backend::template run<scorers>...};
  }
} // namespace mudock
