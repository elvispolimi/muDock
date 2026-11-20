#pragma once

#include <atomic>
#include <chrono>
#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mudock::observability {

  /**
   * Point-in-time snapshot of the live run metrics (from backend_metrics) for a single backend.
   * All fields are plain values (no atomics), so instances are freely copyable and diffable by the observer.
   */
  struct backend_metrics_snapshot {
    std::string_view name;
    std::size_t num_screeners{0};
    std::size_t ligands_batched{0};
    std::size_t batches_emitted{0};
    std::size_t ligands_screening_in_flight{0};
    std::size_t batches_screened{0};
    std::size_t ligands_screened{0};
  };

  /**
   * Point-in-time snapshot of the live run metrics (from run_metrics).
   * All fields are plain values (no atomics), so instances are freely copyable and diffable by the observer.
   */
  struct run_metrics_snapshot {
    std::chrono::steady_clock::time_point time{};
    std::chrono::duration<double> elapsed{};

    std::size_t in_flight{0};
    std::size_t peak_in_flight{0};

    std::size_t ligands_read{0};
    std::size_t ligands_parsed{0};
    std::size_t ligands_dropped_on_parse_fail{0};
    std::vector<backend_metrics_snapshot> backends;
    std::size_t ligands_written{0};

    [[nodiscard]] std::chrono::duration<double> elapsed_since(const run_metrics_snapshot& earlier) const {
      return {time - earlier.time};
    }
  };

  /**
   * Monotonic atomic counter, with relaxed memory ordering.
   */
  class counter {
    std::atomic<std::size_t> value_{0};

  public:
    void increase(std::size_t n = 1) { value_.fetch_add(n, std::memory_order_relaxed); }

    [[nodiscard]] std::size_t load() const { return value_.load(std::memory_order_relaxed); }
  };

  /**
   * Non-monotonic gauge paired with a high-water mark, with relaxed memory ordering.
   */
  class gauge {
    std::atomic<std::size_t> value_{0};
    std::atomic<std::size_t> peak_{0};

  public:
    void increase(std::size_t n = 1) {
      const auto curr      = value_.fetch_add(n, std::memory_order_relaxed) + n;
      std::size_t old_peak = peak_.load(std::memory_order_relaxed);
      while (curr > old_peak && !peak_.compare_exchange_weak(old_peak, curr, std::memory_order_relaxed)) {}
    }

    void decrease(std::size_t n = 1) { value_.fetch_sub(n, std::memory_order_relaxed); }

    [[nodiscard]] std::size_t load() const { return value_.load(std::memory_order_relaxed); }
    [[nodiscard]] std::size_t peak() const { return peak_.load(std::memory_order_relaxed); }
  };

  /**
   * Per-backend, whole-run metrics slice.
   */
  struct backend_metrics {
    std::string name;
    std::size_t num_screeners;
    counter ligands_batched;
    counter batches_emitted;
    gauge ligands_screening_in_flight;
    counter batches_screened;
    counter ligands_screened;

    backend_metrics(std::string n, std::size_t screeners): name(std::move(n)), num_screeners(screeners) {}
  };

  /**
   * Whole-run metrics aggregate with atomic counters for concurrent updates.
   */
  struct run_metrics {
    std::chrono::steady_clock::time_point start_time{std::chrono::steady_clock::now()};

    gauge in_flight;

    counter ligands_read;
    counter ligands_parsed;
    counter ligands_dropped_on_parse_fail;
    std::vector<std::unique_ptr<backend_metrics>> backends;
    counter ligands_written;

    /// Fill a point-in-time snapshot of the live metrics in place.
    void snapshot_into(run_metrics_snapshot& out) const {
      const auto now                    = std::chrono::steady_clock::now();
      out.time                          = now;
      out.elapsed                       = now - start_time;
      out.in_flight                     = in_flight.load();
      out.peak_in_flight                = in_flight.peak();
      out.ligands_read                  = ligands_read.load();
      out.ligands_parsed                = ligands_parsed.load();
      out.ligands_dropped_on_parse_fail = ligands_dropped_on_parse_fail.load();
      out.ligands_written               = ligands_written.load();

      out.backends.resize(backends.size());
      for (std::size_t i = 0; i < backends.size(); ++i) {
        const auto& bm                              = *backends[i];
        out.backends[i].name                        = std::string_view{bm.name};
        out.backends[i].num_screeners               = bm.num_screeners;
        out.backends[i].ligands_batched             = bm.ligands_batched.load();
        out.backends[i].batches_emitted             = bm.batches_emitted.load();
        out.backends[i].ligands_screening_in_flight = bm.ligands_screening_in_flight.load();
        out.backends[i].batches_screened            = bm.batches_screened.load();
        out.backends[i].ligands_screened            = bm.ligands_screened.load();
      }
    }

    /// Capture a point-in-time snapshot of the live metrics.
    [[nodiscard]] run_metrics_snapshot snapshot() const {
      run_metrics_snapshot s;
      snapshot_into(s);
      return s;
    }
  };

} // namespace mudock::observability
