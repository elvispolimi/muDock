#pragma once

#include <atomic>
#include <istream>
#include <limits>
#include <mudock/format/supported_format.hpp>
#include <string>
#include <oneapi/tbb/parallel_pipeline.h>

namespace mudock {

  template<supported_format format>
  class stream_filter {
    std::istream& stream_;
    std::size_t end_;
    const std::size_t max_bytes_per_token_;
    mutable std::string buffered_text_;
    mutable type_of_format<format> format_splitter_;
    mutable bool flushed_ = false;
    std::atomic<bool>* stop_requested;

  public:
    explicit stream_filter(std::istream& in,
                          std::size_t max_bytes_per_token,
                          std::size_t end = std::numeric_limits<std::size_t>::max(),
                          std::atomic<bool>* stop = nullptr);

    // Source filter: reads from the input stream and produces strings
    std::string operator()(oneapi::tbb::flow_control& fc) const;
  };

} // namespace mudock 
