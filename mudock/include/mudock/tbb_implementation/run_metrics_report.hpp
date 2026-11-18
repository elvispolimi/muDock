#pragma once

#include <chrono>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/tbb_implementation/run_metrics.hpp>
#include <ostream>
#include <sstream>
#include <string>

namespace mudock::observability {

  namespace detail {

    /// RAII guard restoring fmtflags/precision on scope exit, so log formatters don't pollute caller stream state.
    class ostream_state_guard {
      std::ostream& out_;
      std::ios::fmtflags flags_;
      std::streamsize prec_;

    public:
      explicit ostream_state_guard(std::ostream& out)
          : out_(out), flags_(out.flags()), prec_(out.precision()) {}
      ~ostream_state_guard() {
        out_.flags(flags_);
        out_.precision(prec_);
      }
      ostream_state_guard(const ostream_state_guard&)            = delete;
      ostream_state_guard& operator=(const ostream_state_guard&) = delete;
    };

  } // namespace detail

  [[nodiscard]] inline double pct(std::size_t n, std::size_t total) {
    return total > 0 ? 100.0 * static_cast<double>(n) / static_cast<double>(total) : 0.0;
  }

  [[nodiscard]] inline double ligands_mb(std::size_t n) {
    constexpr double bytes_per_mb = 1024.0 * 1024.0;
    return static_cast<double>(n * sizeof(static_molecule)) / bytes_per_mb;
  }

  /// Runtime snapshot: instantaneous distribution of in-flight ligands across pipeline buffers.
  inline void log_runtime(std::ostream& out, const run_metrics_snapshot& s) {
    detail::ostream_state_guard guard(out);
    out << std::fixed << std::setprecision(2);

    const auto elapsed = s.elapsed.count();

    std::size_t total_batched  = 0;
    std::size_t total_screened = 0;
    for (const auto& b: s.backends) {
      total_batched += b.ligands_batched;
      total_screened += b.ligands_screened;
    }

    // Derived backlogs
    const auto saturating_sub = [](std::size_t a, std::size_t b) { return a > b ? a - b : std::size_t{0}; };
    const auto parser_in = saturating_sub(s.ligands_read, s.ligands_parsed + s.ligands_dropped_on_parse_fail);
    const auto parsed_buffer = saturating_sub(s.ligands_parsed, total_batched);
    const auto writer_in     = saturating_sub(total_screened, s.ligands_written);

    out << "[observer :: runtime t=" << elapsed << "s]"
        << "  in-flight=" << s.in_flight << " (~" << ligands_mb(s.in_flight) << " MB)"
        << "  peak_in_flight=" << s.peak_in_flight << " (~" << ligands_mb(s.peak_in_flight) << " MB)  |"
        << "  parser_in=" << parser_in << "  parsed_buffer=" << parsed_buffer << "  writer_in=" << writer_in
        << "\n";

    for (const auto& b: s.backends) {
      const auto backend_in_flight = saturating_sub(b.ligands_batched, b.ligands_screened);
      const auto screening         = b.ligands_screening_in_flight;
      const auto waiting           = saturating_sub(backend_in_flight, screening);
      out << "  [" << b.name << "]  in-flight=" << backend_in_flight << " (~" << ligands_mb(backend_in_flight)
          << " MB)  waiting=" << waiting << "  screening=" << screening << "\n";
    }
  }

  /// Convenience overload to print runtime snapshot via `info()`.
  inline void log_runtime(const run_metrics_snapshot& s) {
    std::ostringstream oss;
    log_runtime(oss, s);
    auto str = oss.str();
    if (!str.empty() && str.back() == '\n') {
      str.pop_back();
    }
    info(str);
  }

  /// Cumulative snapshot: counters since run start, each with the whole-run average throughput.
  inline void log_cumulative(std::ostream& out, const run_metrics_snapshot& s) {
    detail::ostream_state_guard guard(out);
    out << std::fixed << std::setprecision(2);

    const auto elapsed     = s.elapsed.count();
    const auto inv_elapsed = elapsed > 0.0 ? 1.0 / elapsed : 0.0;
    const auto avg         = [inv_elapsed](std::size_t n) { return static_cast<double>(n) * inv_elapsed; };

    const auto parse_failed_pct = pct(s.ligands_dropped_on_parse_fail, s.ligands_read);

    std::size_t total_screened = 0;
    for (const auto& b: s.backends) total_screened += b.ligands_screened;
    const auto inv_screened_pct = total_screened > 0 ? 100.0 / static_cast<double>(total_screened) : 0.0;

    out << "[observer :: totals t=" << elapsed << "s]"
        << "  read=" << s.ligands_read << " (" << avg(s.ligands_read) << "/s)"
        << "  parsed=" << s.ligands_parsed << " (" << avg(s.ligands_parsed) << "/s)"
        << "  dropped_on_parse_fail=" << s.ligands_dropped_on_parse_fail << " (" << parse_failed_pct << "%)"
        << "  written=" << s.ligands_written << " (" << avg(s.ligands_written) << "/s)\n";

    for (const auto& b: s.backends) {
      const auto share_pct     = static_cast<double>(b.ligands_screened) * inv_screened_pct;
      const auto avg_batch     = b.batches_screened > 0 ? static_cast<double>(b.ligands_screened) /
                                                          static_cast<double>(b.batches_screened)
                                                        : 0.0;
      const auto screened_rate = avg(b.ligands_screened);
      const auto per_screener_rate =
          b.num_screeners > 0 ? screened_rate / static_cast<double>(b.num_screeners) : 0.0;
      out << "  [" << b.name << "]"
          << "  share=" << share_pct << "%"
          << "  batched=" << b.ligands_batched << " (" << avg(b.ligands_batched) << "/s)"
          << "  batches_emitted=" << b.batches_emitted << " (" << avg(b.batches_emitted) << "/s)"
          << "  batches_screened=" << b.batches_screened << " (" << avg(b.batches_screened) << "/s)"
          << "  ligands_screened=" << b.ligands_screened << " (" << screened_rate << "/s, "
          << per_screener_rate << "/s/screener)"
          << "  avg_batch=" << avg_batch << "\n";
    }
  }

  /// Convenience overload to print cumulative snapshot via `info()`.
  inline void log_cumulative(const run_metrics_snapshot& s) {
    std::ostringstream oss;
    log_cumulative(oss, s);
    auto str = oss.str();
    if (!str.empty() && str.back() == '\n') {
      str.pop_back();
    }
    info(str);
  }

  /// Rolling snapshot: deltas since last snapshot, each with the interval throughput.
  inline void
      log_rolling(std::ostream& out, const run_metrics_snapshot& curr, const run_metrics_snapshot& prev) {
    detail::ostream_state_guard guard(out);
    out << std::fixed << std::setprecision(2);

    const auto dt_sec = curr.elapsed_since(prev).count();
    const auto inv_dt = dt_sec > 0.0 ? 1.0 / dt_sec : 0.0;
    const auto rate   = [inv_dt](std::size_t n) { return static_cast<double>(n) * inv_dt; };

    const auto d_read           = curr.ligands_read - prev.ligands_read;
    const auto d_parsed         = curr.ligands_parsed - prev.ligands_parsed;
    const auto d_parse_failed   = curr.ligands_dropped_on_parse_fail - prev.ligands_dropped_on_parse_fail;
    const auto parse_failed_pct = pct(d_parse_failed, d_read);
    const auto d_written        = curr.ligands_written - prev.ligands_written;

    out << "[observer :: snap Δ=" << dt_sec << "s]"
        << "  read=" << d_read << " (" << rate(d_read) << "/s)"
        << "  parsed=" << d_parsed << " (" << rate(d_parsed) << "/s)"
        << "  dropped_on_parse_fail=" << d_parse_failed << " (" << parse_failed_pct << "%)"
        << "  written=" << d_written << " (" << rate(d_written) << "/s)\n";

    for (std::size_t i = 0; i < curr.backends.size(); ++i) {
      const auto& cb           = curr.backends[i];
      const auto& pb           = prev.backends[i];
      const auto d_batched     = cb.ligands_batched - pb.ligands_batched;
      const auto d_bemitted    = cb.batches_emitted - pb.batches_emitted;
      const auto d_bscreened   = cb.batches_screened - pb.batches_screened;
      const auto d_screened    = cb.ligands_screened - pb.ligands_screened;
      const auto screened_rate = rate(d_screened);
      const auto per_screener_rate =
          cb.num_screeners > 0 ? screened_rate / static_cast<double>(cb.num_screeners) : 0.0;
      out << "  [" << cb.name << "]"
          << "  batched=" << d_batched << " (" << rate(d_batched) << "/s)"
          << "  batches_emitted=" << d_bemitted << " (" << rate(d_bemitted) << "/s)"
          << "  batches_screened=" << d_bscreened << " (" << rate(d_bscreened) << "/s)"
          << "  ligands_screened=" << d_screened << " (" << screened_rate << "/s, " << per_screener_rate
          << "/s/screener)\n";
    }
  }

  /// Convenience overload to print rolling snapshot via `info()`.
  inline void log_rolling(const run_metrics_snapshot& curr, const run_metrics_snapshot& prev) {
    std::ostringstream oss;
    log_rolling(oss, curr, prev);
    auto s = oss.str();
    if (!s.empty() && s.back() == '\n')
      s.pop_back();
    info(s);
  }

  /// Cumulative end-of-run summary: top-level counters + averaged throughput + per-backend breakdown.
  inline void log_summary(std::ostream& out, const run_metrics_snapshot& s) {
    detail::ostream_state_guard guard(out);
    out << std::fixed << std::setprecision(2);

    const auto elapsed     = s.elapsed.count();
    const auto inv_elapsed = elapsed > 0.0 ? 1.0 / elapsed : 0.0;
    const auto avg         = [inv_elapsed](std::size_t n) { return static_cast<double>(n) * inv_elapsed; };

    constexpr int label_w = 22;
    constexpr int value_w = 10;
    const auto count_row  = [&](std::string_view label, std::size_t count) -> std::ostream& {
      out << "  " << std::setw(label_w) << std::left << label << std::setw(value_w) << std::right << count;
      return out;
    };

    out << "Run summary\n"
        << "  " << std::setw(label_w) << std::left << "elapsed" << std::setw(value_w) << std::right << elapsed
        << " s\n";
    count_row("read", s.ligands_read) << "  (" << avg(s.ligands_read) << "/s)\n";
    count_row("parsed", s.ligands_parsed) << "  (" << avg(s.ligands_parsed) << "/s)\n";
    count_row("dropped_on_parse_fail", s.ligands_dropped_on_parse_fail)
        << "  (" << pct(s.ligands_dropped_on_parse_fail, s.ligands_read) << "%)\n";

    if (!s.backends.empty()) {
      std::size_t total_screened = 0;
      for (const auto& b: s.backends) total_screened += b.ligands_screened;
      const auto inv_screened_pct = total_screened > 0 ? 100.0 / static_cast<double>(total_screened) : 0.0;

      // Inner column widths chosen so values right-align at the same column as the top-level ones.
      //   top-level: indent(2) + label_w(22) + value_w(10) == 34
      //   backend:   prefix(9) + label_w(15) + value_w(10) == 34
      constexpr int b_label_w = 15;
      constexpr int b_value_w = 10;

      out << "  backends:\n";
      for (std::size_t i = 0; i < s.backends.size(); ++i) {
        const auto& b        = s.backends[i];
        const bool is_last   = (i + 1 == s.backends.size());
        const char* prefix   = is_last ? "         " : "    │    ";
        const auto avg_batch = b.batches_screened > 0 ? static_cast<double>(b.ligands_screened) /
                                                            static_cast<double>(b.batches_screened)
                                                      : 0.0;
        const auto share_pct = static_cast<double>(b.ligands_screened) * inv_screened_pct;

        const auto screener_throughput =
            b.num_screeners > 0 ? avg(b.ligands_screened) / static_cast<double>(b.num_screeners) : 0.0;

        out << "    " << (is_last ? "└─" : "├─") << " [" << b.name << "]  (" << b.num_screeners
            << " screeners)\n"
            << prefix << std::setw(b_label_w) << std::left << "share" << std::setw(b_value_w) << std::right
            << share_pct << "%\n"
            << prefix << std::setw(b_label_w) << std::left << "batches_screened" << std::setw(b_value_w)
            << std::right << b.batches_screened << "  (" << avg(b.batches_screened) << "/s)\n"
            << prefix << std::setw(b_label_w) << std::left << "avg_batch" << std::setw(b_value_w)
            << std::right << avg_batch << "\n"
            << prefix << std::setw(b_label_w) << std::left << "ligands_screened" << std::setw(b_value_w)
            << std::right << b.ligands_screened << "  (" << avg(b.ligands_screened) << "/s, "
            << screener_throughput << "/s/screener)\n";
      }
    }

    count_row("written", s.ligands_written) << "  (" << avg(s.ligands_written) << "/s)\n";
    out << "  " << std::setw(label_w) << std::left << "peak_in_flight" << std::setw(value_w) << std::right
        << s.peak_in_flight << "  (~" << ligands_mb(s.peak_in_flight) << " MB)\n";
  }

  /// Convenience overload: print summary via `info()`.
  inline void log_summary(const run_metrics_snapshot& s) {
    std::ostringstream oss;
    log_summary(oss, s);
    auto str = oss.str();
    if (!str.empty() && str.back() == '\n') {
      str.pop_back();
    }
    info(str);
  }

} // namespace mudock::observability
