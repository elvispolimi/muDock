#pragma once

#include <cstddef>
#include <istream>
#include <limits>
#include <memory>
#include <mudock/format/supported_format.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <mudock/tbb_implementation/batcher_body.hpp>
#include <mudock/tbb_implementation/parser_body.hpp>
#include <mudock/tbb_implementation/reader_body.hpp>
#include <mudock/tbb_implementation/runtime_services.hpp>
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

    reader_in_t ligand_reader(g, reader_body<format>{knobs.max_bytes_per_token, in, end});
    // todo - ^ handle stop_signal

    global_ligands_limiter_t global_ligands_limiter(g, knobs.max_ligands_in_flight);

    parser_mfn_t ligand_parser(g, tbb_flow::unlimited, parser_body<format>{global_ligands_limiter});
    // todo - ^ handle stop_signal

    parsed_ligands_buf_t parsed_ligands_buf(g);

    writer_state writer_state;
    writer_fn_t writer(g, tbb_flow::serial, writer_body{writer_state, global_ligands_limiter});

    std::vector<device_subgraph> devices;
    devices.reserve(configurations.size());

    for (const auto& configuration: configurations) {
      // todo - properly map "configurations" into subgraph configuration

      auto state = std::make_unique<batcher_state>();

      auto batcher = std::make_unique<batcher_mfn_t>(g, tbb_flow::serial, batcher_body{*state});

      auto batch_buf = std::make_unique<batch_buf_t>(g);

      auto scorer = std::make_unique<scorer_mfn_t>(g, /*concurrency=*/1, scorer_body{});

      tbb_flow::make_edge(parsed_ligands_buf, *batcher);
      tbb_flow::make_edge(tbb_flow::output_port<0>(*batcher), *batch_buf);
      tbb_flow::make_edge(*batch_buf, *scorer);
      tbb_flow::make_edge(tbb_flow::output_port<0>(*scorer), writer);

      devices.push_back({std::move(state), std::move(batcher), std::move(batch_buf), std::move(scorer)});
    }

    tbb_flow::make_edge(ligand_reader, global_ligands_limiter);
    tbb_flow::make_edge(global_ligands_limiter, ligand_parser);
    tbb_flow::make_edge(tbb_flow::output_port<0>(ligand_parser), parsed_ligands_buf);

    // --- Phase 1: drain raw ligands
    ligand_reader.activate();
    g.wait_for_all();

    // --- Phase 2: flush rob partial batches
    for (auto& device: devices) {
      device.state->flush_all([&](lig_batch_t&& batch) {
        device.batch_buf->try_put(std::make_shared<lig_batch_t>(std::move(batch)));
      });
    }
    g.wait_for_all();

    // --- Phase 3: flush writer buffer
    writer_state.flush();

    // todo - observer
    // todo - timeout
  }

} // namespace mudock
