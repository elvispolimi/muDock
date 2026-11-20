#pragma once

#include "stream_engine_types.hpp"

#include <exception>
#include <memory>
#include <mudock/format/reader.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/molecule.hpp>
#include <mudock/tbb_implementation/run_metrics.hpp>
#include <oneapi/tbb/flow_graph.h>
#include <string>

namespace mudock::stream_engine {

  /**
   * Multifunction-node body that parses one input token into a `static_molecule`;
   * emits on port 0, or decrements the global limiter if parsing is skipped/aborted.
   */
  template<supported_format format>
  class parser_body {
    global_ligands_ln& global_ligands_limiter_;
    observability::run_metrics& metrics_;

  public:
    explicit parser_body(global_ligands_ln& global_ligands_limiter, observability::run_metrics& metrics)
        : global_ligands_limiter_(global_ligands_limiter), metrics_(metrics) {}

    void operator()(const std::string& token, parser_mfn::output_ports_type& out) const {
      metrics_.in_flight.increase();
      try {
        auto ligand = std::make_shared<static_molecule>(mudock::parser<format, static_molecule>(token));
        metrics_.ligands_parsed.increase();
        std::get<0>(out).try_put(std::move(ligand));
      } catch (const std::exception&) {
        metrics_.ligands_dropped_on_parse_fail.increase();
        metrics_.in_flight.decrease();
        global_ligands_limiter_.decrementer().try_put(tbb_flow::continue_msg{});
      }
    }
  };

} // namespace mudock::stream_engine
