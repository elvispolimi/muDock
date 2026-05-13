#include <mudock/tbb_implementation/stream_filter.hpp>

namespace mudock {

  template<supported_format format>
  stream_filter<format>::stream_filter(std::istream& in,
                                       std::size_t max_bytes_per_token,
                                       std::size_t end,
                                       std::atomic<bool>* stop)
      : stream_(in),
        max_bytes_per_token_(max_bytes_per_token),
        stop_requested(stop) {
    if (end != std::numeric_limits<std::size_t>::max()) {
      end_ = end;
      return;
    }

    // default case (no MPI involved)
    const std::streampos cur = stream_.tellg();
    stream_.seekg(0, std::ios::end);
    end_ = static_cast<std::size_t>(stream_.tellg());
    stream_.seekg(cur);
  }

  template<supported_format format>
  std::string stream_filter<format>::operator()(oneapi::tbb::flow_control& fc) const {
    auto should_stop = [&]() {
      return stop_requested != nullptr && stop_requested->load(std::memory_order_relaxed);
    };

    while (true) {
      if (should_stop()) {
        buffered_text_.clear();
        flushed_ = true;
        fc.stop();
        return {};
      }

      if (!buffered_text_.empty()) {
        auto input_view = std::string_view{buffered_text_};
        const auto next = format_splitter_.next_molecule_start_index(input_view);
        if (next != std::string_view::npos) {
          std::string token = buffered_text_.substr(0, next);
          buffered_text_.erase(0, next);
          return token;
        }
      }

      const std::streampos pos = stream_.tellg();
      if (pos == std::streampos(-1)) {
        if (!flushed_ && !buffered_text_.empty()) {
          flushed_ = true;
          auto token = std::move(buffered_text_);
          buffered_text_.clear();
          return token;
        }
        fc.stop();
        return {};
      }

      const std::size_t cur = static_cast<std::size_t>(pos);
      if (cur >= end_) {
        if (!flushed_ && !buffered_text_.empty()) {
          flushed_ = true;
          auto token = std::move(buffered_text_);
          buffered_text_.clear();
          return token;
        }
        fc.stop();
        return {};
      }

      const std::size_t bytes_per_token = std::min(max_bytes_per_token_, end_ - cur);
      std::string buf(bytes_per_token, '\0');
      stream_.read(buf.data(), static_cast<std::streamsize>(buf.size()));
      if (stream_.gcount() <= 0) {
        if (!flushed_ && !buffered_text_.empty()) {
          flushed_ = true;
          auto token = std::move(buffered_text_);
          buffered_text_.clear();
          return token;
        }
        fc.stop();
        return {};
      }

      buf.resize(static_cast<std::size_t>(stream_.gcount()));
      buffered_text_ += buf;
    }
  }

  template class stream_filter<supported_format::ADTMOL2>;
  template class stream_filter<supported_format::MOL2>;

} // namespace mudock
