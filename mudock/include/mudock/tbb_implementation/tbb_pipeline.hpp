#pragma once

#include <atomic>
#include <istream>
#include <limits>
#include <mudock/format/supported_format.hpp>
#include <mudock/mudock.hpp>
#include <optional>
#include <string>
#include <vector>

namespace mudock {

  template<supported_format format, typename pipeline_t>
  void run_tbb_pipeline(std::istream& in,
                        const std::vector<std::string>& configurations,
                        const knobs& knobs,
                        pipeline_t& pipeline,
                        std::size_t end = std::numeric_limits<std::size_t>::max(),
                        std::optional<double> time_limit_sec = std::nullopt,
                        std::optional<double> observer_sec = std::nullopt);

} // namespace mudock
