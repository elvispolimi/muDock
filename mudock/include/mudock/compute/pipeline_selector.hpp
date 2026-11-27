#pragma once

#include "mudock/format/supported_format.hpp"

#include <algorithm>
#include <array>
#include <mudock/compute/algorithm.hpp>
#include <mudock/compute/pipeline_registry.hpp>
#include <stdexcept>

namespace mudock {

  template<scoring_function scoring>
  struct scoring_function_format;

  template<typename Derived>
  struct scoring_function_format_base {
    static constexpr bool supports(supported_format input) {
      auto formats = Derived::formats;

      return std::find(formats.begin(), formats.end(), input) != formats.end();
    }
  };

  template<>
  struct scoring_function_format<scoring_function::ADT>
      : scoring_function_format_base<scoring_function_format<scoring_function::ADT>> {
    inline static constexpr std::array formats{supported_format::ADTMOL2, supported_format::PDBQT, supported_format::MOL2};
  };

  template<template<scoring_function> class pipeline_type, scoring_function... scorers, typename callback_type>
  bool dispatch_pipeline_registry(pipeline_registry<pipeline_type<scorers>...>,
                                 const scoring_function requested_scoring,
                                 const supported_format format,
                                 callback_type& callback) {
    bool dispatched = false;
    auto try_scoring = [&]<scoring_function selected_scoring>() {
      if (dispatched || requested_scoring != selected_scoring)
        return;
      if (!scoring_function_format<selected_scoring>::supports(format))
        throw std::runtime_error("File format not supported by the selected scoring function");
      callback.template operator()<pipeline_type<selected_scoring>>();
      dispatched = true;
    };
    (try_scoring.template operator()<scorers>(), ...);
    return dispatched;
  }

  template<typename callback_t>
  void dispatch_selected_pipeline(const search_algorithm search,
                                  const scoring_function scoring,
                                  const supported_format format,
                                  callback_t&& callback) {
    bool dispatched = false;
    if (search == search_algorithm::NONE)
      dispatched = dispatch_pipeline_registry(scoring_pipeline_registry{}, scoring, format, callback);
    else if (search == search_algorithm::GENETIC)
      dispatched = dispatch_pipeline_registry(genetic_pipeline_registry{}, scoring, format, callback);

    if (!dispatched)
      throw std::runtime_error("Selected search/scoring pipeline is not registered");
  }
} // namespace mudock
