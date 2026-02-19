#pragma once

#include <istream>
#include <string>
#include <oneapi/tbb/parallel_pipeline.h>

namespace mudock {

  static constexpr std::size_t max_bytes_per_slice = 1048576; // 1MB

  class stream_filter {
    std::istream& stream_;
    size_t end_;

  public:
    explicit stream_filter(std::istream& in, std::size_t end = std::numeric_limits<std::size_t>::max());
    
    // Source filter: reads from the input stream and produces strings
    std::string operator()(oneapi::tbb::flow_control& fc, 
                           std::size_t max_bytes = max_bytes_per_slice) const;
  };

} // namespace mudock 
