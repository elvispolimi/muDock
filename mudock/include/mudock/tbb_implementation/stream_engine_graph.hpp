#pragma once

#include <atomic>
#include <cstddef>
#include <limits>
#include <memory>
#include <mudock/format/supported_format.hpp>
#include <mudock/knobs.hpp>
#include <mudock/log.hpp>
#include <mudock/tbb_implementation/backend_subgraph.hpp>
#include <mudock/tbb_implementation/parser_body.hpp>
#include <mudock/tbb_implementation/reader_body.hpp>
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

    std::atomic<std::size_t> skipped_ligands{0};

    // --- [SETUP] - Create graph nodes
    reader_in ligand_reader(g, reader_body<format>{knobs.max_bytes_per_token, in, end});

    global_ligands_ln global_ligands_limiter(g, knobs.max_ligands_in_flight);

    parser_mfn ligand_parser(g,
                             tbb_flow::unlimited,
                             parser_body<format>{global_ligands_limiter, skipped_ligands});

    parsed_ligands_buf parsed_ligands_buffer(g);

    auto backends = build_backend_subgraphs(g, knobs, pipeline, configurations);

    writer_state writer_state;
    writer_fn writer(g, tbb_flow::serial, writer_body{writer_state, global_ligands_limiter});

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
    writer_state.flush();

    // todo - observer
  }

} // namespace mudock
