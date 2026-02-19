#pragma once

#include <istream>
#include <limits>
#include <mudock/mudock.hpp>
#include <string>
#include <vector>

namespace mudock {

    static constexpr std::size_t max_tbb_tokens = 4;
    
    void run_tbb_pipeline(std::istream& in,
             const std::vector<std::string>& configurations,
             const knobs& knobs,
             genetic_adt_pipeline& pipeline, 
             std::size_t end = std::numeric_limits<std::size_t>::max());
    
} // namespace mudock
