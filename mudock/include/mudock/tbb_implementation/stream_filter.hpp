#pragma once

#include <istream>
#include <string>
#include <oneapi/tbb/parallel_pipeline.h>

namespace mudock {

  class stream_filter {
    std::istream& stream_;
    size_t end_;
    const size_t max_bytes_per_token_;

  public:
    explicit stream_filter(std::istream& in, 
                          std::size_t max_bytes_per_token, 
                          std::size_t end = std::numeric_limits<std::size_t>::max());
    
    // Source filter: reads from the input stream and produces strings
    std::string operator()(oneapi::tbb::flow_control& fc) const;
  };

} // namespace mudock 
