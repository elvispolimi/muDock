#pragma once

#include "stream_engine_types.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <istream>
#include <limits>
#include <mudock/format/supported_format.hpp>
#include <optional>
#include <string>
#include <string_view>

namespace mudock::stream_engine {

  /**
   * Input-node body that reads bytes from the stream, splits on format-specific molecule boundaries,
   * and emits one token at a time; flushes the remainder on EOF, end-position, or stop.
   */
  template<supported_format format>
  class reader_body {
    const std::size_t max_bytes_per_token_;
    std::istream& stream_;
    std::size_t end_;
    std::atomic<bool>& stop_requested_;

    // Reader state: mutable because of the TBB input-node `operator() const`.
    mutable std::string buffered_text_;
    mutable type_of_format<format> format_splitter_;
    mutable bool flushed_ = false;

  public:
    explicit reader_body(std::size_t max_bytes_per_token,
                         std::istream& in,
                         std::atomic<bool>& stop_requested,
                         std::size_t end = std::numeric_limits<std::size_t>::max())
        : max_bytes_per_token_(max_bytes_per_token), stream_(in), stop_requested_(stop_requested) {
      if (end != std::numeric_limits<std::size_t>::max()) {
        end_ = end;
        return;
      }

      const std::streampos cur = stream_.tellg();
      stream_.seekg(0, std::ios::end);
      end_ = static_cast<std::size_t>(stream_.tellg());
      stream_.seekg(cur);
    }

    std::string operator()(oneapi::tbb::flow_control& fc) const {
      while (true) {
        if (stop_requested_.load(std::memory_order_relaxed)) {
          buffered_text_.clear();
          flushed_ = true;
          fc.stop();
          return {};
        }
        if (auto token = try_emit_token()) {
          return std::move(*token);
        }
        if (!read_more()) {
          return flush_remainder(fc);
        }
      }
    }

  private:
    std::optional<std::string> try_emit_token() const {
      if (!buffered_text_.empty()) {
        const auto next_token_pos = format_splitter_.next_molecule_start_index(buffered_text_);
        if (next_token_pos != std::string_view::npos) {
          std::string token = buffered_text_.substr(0, next_token_pos);
          buffered_text_.erase(0, next_token_pos);
          return token;
        }
      }
      return std::nullopt;
    }

    bool read_more() const {
      // Position-tracking unavailable: failed state.
      const std::streampos pos = stream_.tellg();
      if (pos == std::streampos(-1)) {
        return false;
      }

      // End position reached: EOF (seekable stream).
      const auto cur = static_cast<std::size_t>(pos);
      if (cur >= end_) {
        return false;
      }

      const auto bytes_per_token = std::min(max_bytes_per_token_, end_ - cur);
      std::string buf(bytes_per_token, '\0');
      stream_.read(buf.data(), static_cast<std::streamsize>(buf.size()));

      // Stream returned 0 bytes: EOF (non-seekable stream), or unexpected failure (seekable stream).
      const auto bytes_read = stream_.gcount();
      if (bytes_read <= 0) {
        return false;
      }

      buffered_text_.append(buf.data(), static_cast<std::size_t>(bytes_read));
      return true;
    }

    std::string flush_remainder(oneapi::tbb::flow_control& fc) const {
      if (!flushed_ && !buffered_text_.empty()) {
        flushed_   = true;
        auto token = std::move(buffered_text_);
        buffered_text_.clear();
        return token;
      }
      fc.stop();
      return {};
    }
  };

} // namespace mudock::stream_engine
