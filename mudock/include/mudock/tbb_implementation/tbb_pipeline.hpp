#pragma once

#include <istream>
#include <limits>
#include <mudock/mudock.hpp>
#include <string>
#include <vector>

namespace mudock {
    
    template<typename pipeline_t>
    void run_tbb_pipeline(std::istream& in,
             const std::vector<std::string>& configurations,
             const knobs& knobs,
             pipeline_t& pipeline, 
             std::size_t end = std::numeric_limits<std::size_t>::max());
    
} // namespace mudock
