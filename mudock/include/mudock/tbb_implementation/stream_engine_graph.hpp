#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <mudock/knobs.hpp>
#include <mudock/log.hpp>
#include <mudock/tbb_implementation/backend_subgraph.hpp>
#include <mudock/tbb_implementation/parser_body.hpp>
#include <mudock/tbb_implementation/reader_body.hpp>
#include <mudock/tbb_implementation/run_metrics.hpp>
#include <mudock/tbb_implementation/run_metrics_report.hpp>
#include <mudock/tbb_implementation/runtime_services.hpp>
#include <mudock/tbb_implementation/stream_engine_types.hpp>
#include <mudock/tbb_implementation/writer_body.hpp>
#include <oneapi/tbb/flow_graph.h>
#include <optional>
#include <string>
#include <vector>

namespace mudock {

  template<supported_format format, typename pipeline_t>
  void run_stream_engine_graph(std::istream& in,
                               const std::vector<std::string>& configurations,
                               const knobs& knobs,
                               pipeline_t& pipeline,
                               std::size_t end                      = std::numeric_limits<std::size_t>::max(),
                               std::optional<double> time_limit_sec = std::nullopt,
                               std::optional<double> observer_sec   = std::nullopt) {
    using namespace stream_engine;

    tbb_flow::graph g;

    observability::run_metrics metrics;

    // --- [SETUP] - Create graph nodes
    reader_in ligand_reader(g, reader_body<format>{knobs.max_bytes_per_token, in, metrics, end});

    global_ligands_ln global_ligands_limiter(g, knobs.max_ligands_in_flight);

    parser_mfn ligand_parser(g, tbb_flow::unlimited, parser_body<format>{global_ligands_limiter, metrics});

    parsed_ligands_buf parsed_ligands_buffer(g);

    auto backends = build_backend_subgraphs(g, knobs, pipeline, configurations, metrics);

    writer_state writing;
    writer_fn writer(g, tbb_flow::serial, writer_body{writing, global_ligands_limiter, metrics});

    // --- [SETUP] - Wire graph nodes
    tbb_flow::make_edge(ligand_reader, global_ligands_limiter);
    tbb_flow::make_edge(global_ligands_limiter, ligand_parser);
    tbb_flow::make_edge(tbb_flow::output_port<0>(ligand_parser), parsed_ligands_buffer);
    for (auto& backend: backends) backend.wire(parsed_ligands_buffer, writer);

    // --- [SETUP] - Start timer for time limit, if specified
    detail::deadline_timer timer;
    if (time_limit_sec && *time_limit_sec > 0.0) {
      info(std::format("Time limit enabled: {} s", *time_limit_sec));
      timer.start(time_limit_sec, [&, limit = *time_limit_sec]() {
        g.cancel();
        info(std::format("Time limit of {} s reached: cancelling graph ...", limit));
      });
    }

    // --- [SETUP] - Start observer for periodic metrics snapshots, if specified
    detail::periodic_observer observer;
    metrics.start_time = std::chrono::steady_clock::now();
    if (observer_sec && *observer_sec > 0.0) {
      info(std::format("Observer enabled: every {} s", *observer_sec));
      auto prev = metrics.snapshot();
      auto curr = observability::run_metrics_snapshot{};
      observer.start(observer_sec, [&metrics, prev, curr]() mutable {
        metrics.snapshot_into(curr);
        observability::log_runtime(curr);
        observability::log_cumulative(curr);
        observability::log_rolling(curr, prev);
        std::swap(prev, curr);
      });
    }

    // --- [SCREENING] - Phase 1: drain raw ligands
    info("Starting virtual screening ...");
    ligand_reader.activate();
    g.wait_for_all();
    timer.cancel();

    // --- [SCREENING] - Phase 2: flush ROBs' partial batches (skipped if time limit expired)
    if (!g.is_cancelled()) {
      info("Flushing partial batches ...");
      for (auto& backend: backends) backend.flush_batches();
      g.wait_for_all();
    }

    // --- [SCREENING] - Phase 3: flush writer's buffer
    info("Flushing writer's buffer ...");
    writing.flush();

    // --- [DONE] - Final run summary
    observer.stop();
    observability::log_summary(metrics.snapshot());
  }

} // namespace mudock
