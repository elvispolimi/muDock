#pragma once

#include <istream>
#include <string>
#include <vector>

#include <mudock/mudock.hpp>

namespace mudock {

    static constexpr std::size_t max_tbb_tokens = 4;
    
    void run_tbb_pipeline(std::istream& in,
             const std::vector<std::string>& configurations,
             const knobs& knobs,
             genetic_adt_pipeline& pipeline, 
             std::size_t max_tokens = max_tbb_tokens);
    
} // namespace mudock
