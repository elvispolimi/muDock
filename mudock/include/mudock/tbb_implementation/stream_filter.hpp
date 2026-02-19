#pragma once

#include <istream>
#include <string>
#include <oneapi/tbb/parallel_pipeline.h>

namespace mudock {

  static constexpr std::size_t max_bytes_per_slice = 20000;

  class stream_filter {
    std::istream& stream_;

  public:
    explicit stream_filter(std::istream& in);
    
    // Source filter: reads from the input stream and produces strings
    std::string operator()(oneapi::tbb::flow_control& fc, 
                           std::size_t max_bytes = max_bytes_per_slice) const;
  };

} // namespace mudock 
