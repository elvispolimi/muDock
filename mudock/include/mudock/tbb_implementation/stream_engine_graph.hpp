#pragma once

#include <atomic>
#include <cstddef>
#include <limits>
#include <memory>
#include <mudock/format/supported_format.hpp>
#include <mudock/knobs.hpp>
#include <mudock/tbb_implementation/batcher_body.hpp>
#include <mudock/tbb_implementation/parser_body.hpp>
#include <mudock/tbb_implementation/reader_body.hpp>
#include <mudock/tbb_implementation/scorer_body.hpp>
#include <mudock/tbb_implementation/stream_engine_types.hpp>
#include <mudock/tbb_implementation/writer_body.hpp>
#include <oneapi/tbb/flow_graph.h>
#include <optional>
#include <string>
#include <utility>
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

    std::atomic<bool> stop_requested{false};
    std::atomic<std::size_t> skipped_ligands{0};

    reader_in ligand_reader(g, reader_body<format>{knobs.max_bytes_per_token, in, stop_requested, end});

    global_ligands_limiter global_ligands_limiter(g, knobs.max_ligands_in_flight);

    parser_mfn ligand_parser(g,
                             tbb_flow::unlimited,
                             parser_body<format>{global_ligands_limiter, skipped_ligands, stop_requested});

    parsed_ligands_buf parsed_ligands_buf(g);

    writer_state writer_state;
    writer_fn writer(g, tbb_flow::serial, writer_body{writer_state, global_ligands_limiter});

    std::vector<device_subgraph> devices;
    devices.reserve(configurations.size());

    for (const auto& configuration: configurations) {
      // todo - properly map "configurations" into subgraph configuration

      auto state = std::make_unique<batcher_state>();

      auto batcher = std::make_unique<batcher_mfn>(g, tbb_flow::serial, batcher_body{*state});

      auto batch_buffer = std::make_unique<batch_buf>(g);

      auto scorer = std::make_unique<scorer_mfn>(g, /*concurrency=*/1, scorer_body{});

      tbb_flow::make_edge(parsed_ligands_buf, *batcher);
      tbb_flow::make_edge(tbb_flow::output_port<0>(*batcher), *batch_buffer);
      tbb_flow::make_edge(*batch_buffer, *scorer);
      tbb_flow::make_edge(tbb_flow::output_port<0>(*scorer), writer);

      devices.push_back({std::move(state), std::move(batcher), std::move(batch_buffer), std::move(scorer)});
    }

    tbb_flow::make_edge(ligand_reader, global_ligands_limiter);
    tbb_flow::make_edge(global_ligands_limiter, ligand_parser);
    tbb_flow::make_edge(tbb_flow::output_port<0>(ligand_parser), parsed_ligands_buf);

    // --- Phase 1: drain raw ligands
    ligand_reader.activate();
    g.wait_for_all();

    // --- Phase 2: flush rob partial batches
    for (auto& device: devices) {
      device.state->flush_all([&](lig_batch&& batch) {
        device.batch_buffer->try_put(std::make_shared<lig_batch>(std::move(batch)));
      });
    }
    g.wait_for_all();

    // --- Phase 3: flush writer buffer
    writer_state.flush();

    // todo - observer
    // todo - timeout
  }

} // namespace mudock
