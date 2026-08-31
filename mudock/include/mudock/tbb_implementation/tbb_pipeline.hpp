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
    std::atomic<bool> measurement_started{false};
    std::atomic<bool> measurement_complete{false};
    std::atomic<bool> measurement_stop_done{false};
    std::atomic<bool> measurement_watcher_exit{false};
    std::atomic<std::size_t> measurement_start_count{0};
    std::atomic<std::size_t> measurement_end_count{0};
    std::atomic<std::size_t> measurement_batches_completed{0};
    std::atomic<std::size_t> measurement_start_batch_count{0};
    std::atomic<std::int64_t> measurement_start_ns{0};
    std::atomic<std::int64_t> measurement_end_ns{0};
    const bool ligand_measurement = measure_ligands && *measure_ligands > 0;
    const bool batch_measurement  = measure_batches && *measure_batches > 0;
    const bool fixed_measurement  = ligand_measurement || batch_measurement;
    const auto start = std::chrono::high_resolution_clock::now();

    const auto timestamp_ns = [] {
      return std::chrono::duration_cast<std::chrono::nanoseconds>(
                 std::chrono::steady_clock::now().time_since_epoch())
          .count();
    };

    const auto complete_if_ready = [&](const std::size_t observed_output_count) {
      if (!measurement_started.load(std::memory_order_relaxed)) {
        return;
      }
      const auto start_count  = measurement_start_count.load(std::memory_order_relaxed);
      const auto end_count    = observed_output_count;
      const auto start_batches = measurement_start_batch_count.load(std::memory_order_relaxed);
      const auto end_batches   = measurement_batches_completed.load(std::memory_order_relaxed);
      const bool ligands_ready = !ligand_measurement ||
                                 (end_count >= start_count && end_count - start_count >= *measure_ligands);
      const bool batches_ready = !batch_measurement ||
                                 (end_batches >= start_batches && end_batches - start_batches >= *measure_batches);
      if (!ligands_ready || !batches_ready) {
        return;
      }

      bool expected = false;
      if (measurement_complete.compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
        measurement_end_count.store(end_count, std::memory_order_relaxed);
        measurement_end_ns.store(timestamp_ns(), std::memory_order_relaxed);
      }
    };

    if (fixed_measurement) {
      info("Steady-state measurement enabled from first submitted batch: measure_ligands=",
           ligand_measurement ? std::to_string(*measure_ligands) : std::string{"disabled"},
           ", measure_batches=",
           batch_measurement ? std::to_string(*measure_batches) : std::string{"disabled"});

      output_queue->set_enqueue_callback([&, complete_if_ready](const std::size_t count) {
        const auto start_count = measurement_start_count.load(std::memory_order_relaxed);
        if (ligand_measurement && measurement_started.load(std::memory_order_relaxed) && count >= start_count &&
            count - start_count >= *measure_ligands) {
          complete_if_ready(count);
        }
      });
    }

    const auto on_batch_submitted = [&] {
      if (!fixed_measurement) {
        return;
      }
      bool expected = false;
      if (measurement_started.compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
        measurement_start_count.store(output_queue->get_global_counter(), std::memory_order_relaxed);
        measurement_start_batch_count.store(measurement_batches_completed.load(std::memory_order_relaxed),
                                             std::memory_order_relaxed);
        measurement_start_ns.store(timestamp_ns(), std::memory_order_relaxed);
      }
    };

    const auto on_batch_completed = [&] {
      const auto completed = measurement_batches_completed.fetch_add(1, std::memory_order_relaxed) + 1;
      const auto start_batches = measurement_start_batch_count.load(std::memory_order_relaxed);
      if (batch_measurement && measurement_started.load(std::memory_order_relaxed) && completed >= start_batches &&
          completed - start_batches >= *measure_batches) {
        complete_if_ready(output_queue->get_global_counter());
      }
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

    // This is deliberately outside safe_queue::enqueue: enqueue callbacks run
    // while the queue mutex is held and therefore cannot close a queue safely.
    auto stop_without_drain = [&]() {
      bool expected = false;
      if (!measurement_stop_done.compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
        return;
      }
      stop_requested.store(true, std::memory_order_relaxed);
      input_queue->send_terminate_signal();
      dropped_by_timeout.store(input_queue->clear(), std::memory_order_relaxed);
      output_queue->clear();
      output_queue->send_terminate_signal();
    };

    std::thread measurement_watcher;
    if (fixed_measurement) {
      measurement_watcher = std::thread([&]() {
        while (!measurement_watcher_exit.load(std::memory_order_relaxed) &&
               !measurement_complete.load(std::memory_order_relaxed) &&
               !timeout_triggered.load(std::memory_order_relaxed)) {
          std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        if (measurement_complete.load(std::memory_order_relaxed) ||
            timeout_triggered.load(std::memory_order_relaxed)) {
          stop_without_drain();
        }
      });
    }

    if (time_limit_sec && *time_limit_sec > 0.0) {
      info("Time limit enabled: ", *time_limit_sec, " s");
      timer.start(time_limit_sec, [&]() {
        stop_requested.store(true, std::memory_order_relaxed);
        input_queue->send_terminate_signal();
        dropped_by_timeout.store(input_queue->clear(), std::memory_order_relaxed);
        timeout_triggered.store(true, std::memory_order_relaxed);
        if (fixed_measurement) {
          stop_without_drain();
        }
        info("Time limit reached: discarded ",
             dropped_by_timeout.load(std::memory_order_relaxed),
             " pending ligand(s) from input queue.");
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
      measurement_watcher_exit.store(true, std::memory_order_relaxed);
      if (measurement_watcher.joinable()) {
        measurement_watcher.join();
      }
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
    if (fixed_measurement && measurement_complete.load(std::memory_order_relaxed)) {
      const auto start_ns = measurement_start_ns.load(std::memory_order_relaxed);
      const auto end_ns   = measurement_end_ns.load(std::memory_order_relaxed);
      const auto count0   = measurement_start_count.load(std::memory_order_relaxed);
      const auto count1   = measurement_end_count.load(std::memory_order_relaxed);
      const double elapsed = static_cast<double>(end_ns - start_ns) * 1.0e-9;
      const double throughput = elapsed > 0.0 ? static_cast<double>(count1 - count0) / elapsed : 0.0;
      info("Steady-state measurement complete: processed=", count1 - count0,
           ", elapsed=", elapsed, " s, throughput=", throughput,
           " ligands/s; stopped without draining pending output.");
    } else if (fixed_measurement) {
      info("Steady-state measurement incomplete: processed=", output_queue->get_global_counter(),
           ", completed_batches=", measurement_batches_completed.load(std::memory_order_relaxed));
    } else {
      info("Output drained");
    }
  }

} // namespace mudock
