#pragma once

#include <atomic>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <istream>
#include <limits>
#include <memory>
#include <mudock/compute/manager.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>
#include <mudock/tbb_implementation/parser_filter.hpp>
#include <mudock/tbb_implementation/runtime_services.hpp>
#include <mudock/tbb_implementation/stream_filter.hpp>
#include <oneapi/tbb/parallel_pipeline.h>
#include <optional>
#include <string>
#include <thread>
#include <vector>

namespace mudock {

  namespace detail {
    inline void append_ligand(std::string& buf, const static_molecule& ligand) {
      buf += ligand.properties.get(property_type::NAME);
      buf += ' ';
      buf += ligand.properties.get(property_type::SCORE);
      buf += '\n';
    }

    struct measurement_state {
      std::atomic<std::int64_t> started_at_ns{-1};
      std::atomic<std::size_t> completed_ligands{0};
      std::atomic<std::size_t> completed_batches{0};
      std::atomic<bool> shutdown_started{false};
    };
  } // namespace detail

  template<supported_format format, typename pipeline_t>
  void run_tbb_pipeline(std::istream& in,
                        const std::vector<std::string>& configurations,
                        const knobs& knobs,
                        pipeline_t& pipeline,
                        std::size_t end = std::numeric_limits<std::size_t>::max(),
                        std::optional<double> time_limit_sec = std::nullopt,
                        std::optional<double> observer_sec = std::nullopt,
                        std::optional<std::size_t> measure_ligands = std::nullopt,
                        std::optional<std::size_t> measure_batches = std::nullopt) {
    auto input_queue  = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
    auto output_queue = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
    input_queue->initialize(knobs.max_tbb_queue_size);
    output_queue->initialize(knobs.max_tbb_queue_size);
    std::atomic<std::size_t> skipped_ligands{0};
    std::atomic<std::size_t> dropped_by_timeout{0};
    std::atomic<std::size_t> in_flight_ligands{0};
    std::atomic<bool> timeout_triggered{false};
    std::atomic<bool> stop_requested{false};
    detail::measurement_state measurement;
    const bool ligand_measurement = measure_ligands && *measure_ligands > 0;
    const bool batch_measurement  = measure_batches && *measure_batches > 0;
    const bool fixed_measurement  = ligand_measurement || batch_measurement;
    const auto start = std::chrono::high_resolution_clock::now();

    const auto timestamp_ns = [] {
      return std::chrono::duration_cast<std::chrono::nanoseconds>(
                 std::chrono::steady_clock::now().time_since_epoch())
          .count();
    };

    // This is deliberately outside safe_queue::enqueue: callbacks run after
    // the queue mutex is released and may therefore close the queues safely.
    const auto stop_without_drain = [&] {
      stop_requested.store(true, std::memory_order_relaxed);
      input_queue->send_terminate_signal();
      dropped_by_timeout.store(input_queue->clear(), std::memory_order_relaxed);
      output_queue->clear();
      output_queue->send_terminate_signal();
    };

    const auto complete_if_ready = [&] {
      const auto completed_ligands = measurement.completed_ligands.load(std::memory_order_relaxed);
      const auto completed_batches = measurement.completed_batches.load(std::memory_order_relaxed);
      const bool ligands_ready = !ligand_measurement || completed_ligands >= *measure_ligands;
      const bool batches_ready = !batch_measurement || completed_batches >= *measure_batches;
      if ((!ligand_measurement && !batch_measurement) || !ligands_ready || !batches_ready)
        return;

      bool expected = false;
      if (!measurement.shutdown_started.compare_exchange_strong(expected, true, std::memory_order_relaxed))
        return;

      const auto end_ns = timestamp_ns();
      const auto start_ns = measurement.started_at_ns.load(std::memory_order_acquire);
      const double elapsed = start_ns >= 0
                                 ? static_cast<double>(end_ns - start_ns) * 1.0e-9
                                 : 0.0;
      const double throughput = elapsed > 0.0 ? static_cast<double>(completed_ligands) / elapsed : 0.0;
      info("Steady-state measurement complete: ligands=",
           completed_ligands,
           ", batches=",
           completed_batches,
           ", elapsed=",
           elapsed,
           " s, throughput=",
           throughput,
           " ligands/s; stopped without draining pending output.");
      stop_without_drain();
    };

    if (fixed_measurement) {
      info("Steady-state measurement enabled from first submitted batch: measure_ligands=",
           ligand_measurement ? std::to_string(*measure_ligands) : std::string{"disabled"},
           ", measure_batches=",
           batch_measurement ? std::to_string(*measure_batches) : std::string{"disabled"});

      output_queue->set_enqueue_callback([&](const std::size_t) {
        if (!ligand_measurement || measurement.started_at_ns.load(std::memory_order_acquire) < 0)
          return;
        measurement.completed_ligands.fetch_add(1, std::memory_order_relaxed);
        complete_if_ready();
      });
    }

    const auto on_batch_submitted = [&] {
      if (!fixed_measurement)
        return;
      std::int64_t expected = -1;
      measurement.started_at_ns.compare_exchange_strong(
          expected, timestamp_ns(), std::memory_order_release, std::memory_order_relaxed);
    };

    const auto on_batch_completed = [&] {
      if (measurement.started_at_ns.load(std::memory_order_acquire) < 0)
        return;
      measurement.completed_batches.fetch_add(1, std::memory_order_relaxed);
      complete_if_ready();
    };

    std::thread writer([&] {
      std::string buf;
      buf.reserve(1 << 20);

      std::size_t lines = 0;
      for (auto x = output_queue->dequeue(); x; x = output_queue->dequeue()) {
        detail::append_ligand(buf, *x);

        if (++lines % 4096 == 0) {
          std::cout << buf;
          buf.clear();
        }
      }
      if (!buf.empty())
        std::cout << buf;
    });

    threadpool pool;

    std::size_t prev_processed = output_queue->get_global_counter();
    auto prev_time             = std::chrono::high_resolution_clock::now();
    detail::periodic_observer observer;
    if (observer_sec && *observer_sec > 0.0) {
      info("Observer enabled with period: ", *observer_sec, " s");
      observer.start(observer_sec, [&]() {
        const auto now                  = std::chrono::high_resolution_clock::now();
        const std::size_t now_processed = output_queue->get_global_counter();
        const std::size_t in_backlog    = input_queue->size();
        const std::size_t in_flight     = in_flight_ligands.load(std::memory_order_relaxed);

        const std::chrono::duration<double> dt = now - prev_time;
        const std::size_t delta_processed      = now_processed - prev_processed;
        const double inst_throughput =
            dt.count() > 0.0 ? static_cast<double>(delta_processed) / dt.count() : 0.0;
        const std::chrono::duration<double> total = now - start;
        const double avg_throughput =
            total.count() > 0.0 ? static_cast<double>(now_processed) / total.count() : 0.0;

        info("Observer: processed=",
             now_processed,
             ", input_backlog=",
             in_backlog,
             ", in_flight=",
             in_flight,
             ", inst_throughput=",
             inst_throughput,
             " ligands/s, avg_throughput=",
             avg_throughput,
             " ligands/s");

        prev_processed = now_processed;
        prev_time      = now;
      });
    }

    detail::deadline_timer timer;

    if (time_limit_sec && *time_limit_sec > 0.0) {
      info("Time limit enabled: ", *time_limit_sec, " s");
      timer.start(time_limit_sec, [&]() {
        bool expected = false;
        if (measurement.shutdown_started.compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
          timeout_triggered.store(true, std::memory_order_relaxed);
          stop_without_drain();
          info("Time limit reached: discarded ",
               dropped_by_timeout.load(std::memory_order_relaxed),
               " pending ligand(s) from input queue.");
        }
      });
    }

    bool finalized = false;
    auto finalize_runtime = [&]() {
      if (finalized) {
        return;
      }
      finalized = true;

      stop_requested.store(true, std::memory_order_relaxed);
      input_queue->send_terminate_signal();
      pool.wait();
      output_queue->send_terminate_signal();

      if (writer.joinable()) {
        writer.join();
      }
      observer.stop();
      observer.join();
      timer.cancel();
      timer.join();
    };

    try {
      manager(configurations,
              pool,
              knobs,
              input_queue,
              output_queue,
              pipeline,
              &in_flight_ligands,
              on_batch_submitted,
              on_batch_completed);
      info("Manager done: workers created");

      oneapi::tbb::parallel_pipeline(
          knobs.max_tbb_tokens,
          oneapi::tbb::make_filter<void, std::string>(
              oneapi::tbb::filter_mode::serial_in_order,
              stream_filter<format>(in, knobs.max_bytes_per_token, end, &stop_requested))
              & oneapi::tbb::make_filter<std::string, void>(
                  // OpenBabel MOL2 parsing is deliberately serialized; the
                  // ADTMOL2/native parser remains parallelizable.
                  format == supported_format::MOL2 ? oneapi::tbb::filter_mode::serial_in_order
                                                   : oneapi::tbb::filter_mode::parallel,
                  parser_filter<format>(input_queue, &skipped_ligands, &stop_requested)));

      info("Pipeline done: closing input queue");
      finalize_runtime();
    } catch (...) {
      finalize_runtime();
      throw;
    }

    if (const auto skipped = skipped_ligands.load(std::memory_order_relaxed); skipped > 0) {
      error("Skipped ", skipped, " ligand(s) due to parse errors.");
    }
    if (timeout_triggered.load(std::memory_order_relaxed)) {
      info("Dropped ligands due to timeout: ", dropped_by_timeout.load(std::memory_order_relaxed));
    }
    const bool measurement_completed = measurement.shutdown_started.load(std::memory_order_relaxed) &&
                                       !timeout_triggered.load(std::memory_order_relaxed);
    if (fixed_measurement && !measurement_completed) {
      const auto start_ns = measurement.started_at_ns.load(std::memory_order_acquire);
      const double elapsed = start_ns >= 0
                                 ? static_cast<double>(timestamp_ns() - start_ns) * 1.0e-9
                                 : 0.0;
      info("Steady-state measurement incomplete: ligands=",
           measurement.completed_ligands.load(std::memory_order_relaxed),
           ", batches=",
           measurement.completed_batches.load(std::memory_order_relaxed),
           ", elapsed=",
           elapsed,
           " s");
    } else if (!fixed_measurement) {
      info("Output drained");
    }
  }

} // namespace mudock
