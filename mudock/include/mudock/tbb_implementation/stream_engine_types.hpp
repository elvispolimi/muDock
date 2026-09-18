#pragma once

#include <memory>
#include <mudock/batch.hpp>
#include <mudock/molecule.hpp>
#include <oneapi/tbb/flow_graph.h>
#include <string>
#include <tuple>

namespace mudock::stream_engine {
  namespace tbb_flow = oneapi::tbb::flow;

  using lig_ptr       = std::shared_ptr<static_molecule>;
  using lig_batch     = batch<static_molecule>;
  using lig_batch_ptr = std::shared_ptr<lig_batch>;

  using reader_in              = tbb_flow::input_node<std::string>;
  using global_ligands_limiter = tbb_flow::limiter_node<std::string>;
  using parser_mfn             = tbb_flow::multifunction_node<std::string, std::tuple<lig_ptr>>;
  using parsed_ligands_buf     = tbb_flow::buffer_node<lig_ptr>;
  using batcher_mfn = tbb_flow::multifunction_node<lig_ptr, std::tuple<lig_batch_ptr>, tbb_flow::rejecting>;
  using batch_buf   = tbb_flow::buffer_node<lig_batch_ptr>;
  using scorer_mfn  = tbb_flow::multifunction_node<lig_batch_ptr, std::tuple<lig_ptr>, tbb_flow::rejecting>;
  using writer_fn   = tbb_flow::function_node<lig_ptr>;

  class batcher_state;

  struct device_subgraph {
    std::unique_ptr<batcher_state> state;
    std::unique_ptr<batcher_mfn> batcher;
    std::unique_ptr<batch_buf> batch_buffer;
    std::unique_ptr<scorer_mfn> scorer;
  };

} // namespace mudock::stream_engine
