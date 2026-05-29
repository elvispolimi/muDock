#pragma once

#include "mudock/format/supported_format.hpp"
#include "mudock/utils.hpp"

#include <algorithm>
#include <array>
#include <mudock/compute/algorithm.hpp>
#include <mudock/compute/pipeline.hpp>
#include <stdexcept>

namespace mudock {

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
    inline static constexpr std::array formats{supported_format::ADTMOL2};
  };

  template<typename callback_t>
  void dispatch_selected_pipeline(const search_algorithm search,
                                  const scoring_function scoring,
                                  const supported_format format,
                                  callback_t&& callback) {
    constexpr_switch<0, num_search_algorithms(), 1>(
        [&](const auto search_index) {
          constexpr auto selected_search = static_cast<search_algorithm>(decltype(search_index)::value);
          constexpr_switch<0, num_scoring_functions(), 1>(
              [&](const auto scoring_index) {
                constexpr auto selected_scoring =
                    static_cast<scoring_function>(decltype(scoring_index)::value);
                if (!scoring_function_format<selected_scoring>::supports(format))
                  throw std::runtime_error("File format not supported by the selected scoring function");
                using pipeline_t = pipeline_selector_t<selected_search, selected_scoring>;
                callback.template operator()<pipeline_t>(selected_search, selected_scoring);
              },
              static_cast<int>(scoring));
        },
        static_cast<int>(search));
  }
} // namespace mudock
