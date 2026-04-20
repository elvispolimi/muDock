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

    auto flush_remainder = [&]() {
      if (flushed_) {
        return;
      }
      flushed_ = true;
      if (!buffered_text_.empty()) {
        ready_tokens_.emplace_back(std::move(buffered_text_));
        buffered_text_.clear();
      }
    };

    while (ready_tokens_.empty()) {
      if (should_stop()) {
        ready_tokens_.clear();
        fc.stop();
        return {};
      }

      const std::streampos pos = stream_.tellg();
      if (pos == std::streampos(-1)) {
        flush_remainder();
        break;
      }

      const std::size_t cur = static_cast<std::size_t>(pos);
      if (cur >= end_) {
        flush_remainder();
        break;
      }

      const std::size_t bytes_per_token = std::min(max_bytes_per_token_, end_ - cur);
      std::string buf(bytes_per_token, '\0');
      stream_.read(buf.data(), static_cast<std::streamsize>(buf.size()));
      if (stream_.gcount() <= 0) {
        flush_remainder();
        break;
      }

      buf.resize(static_cast<std::size_t>(stream_.gcount()));
      buffered_text_ += buf;
      auto input_view = std::string_view{buffered_text_};
      while (!input_view.empty()) {
        const auto next = format_splitter_.next_molecule_start_index(input_view);
        if (next == std::string_view::npos) {
          buffered_text_ = std::string{input_view};
          break;
        }

        ready_tokens_.emplace_back(input_view.substr(0, next));
        input_view = input_view.substr(next);
      }
    }

    if (ready_tokens_.empty()) {
      fc.stop();
      return {};
    }

    auto token = std::move(ready_tokens_.front());
    ready_tokens_.pop_front();
    return token;
  }

  template class stream_filter<supported_format::ADTMOL2>;
  template class stream_filter<supported_format::MOL2>;

} // namespace mudock
