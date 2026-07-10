#pragma once

#include "stream_engine_types.hpp"

#include <cstddef>
#include <iostream>
#include <mudock/molecule.hpp>
#include <mudock/molecule/properties.hpp>
#include <oneapi/tbb/flow_graph.h>
#include <string>

namespace mudock::stream_engine {

  class writer_state {
    static constexpr std::size_t flush_every    = 4096;
    static constexpr std::size_t buffer_reserve = 1 << 20;

    std::string buffer_;
    std::size_t lines_ = 0;

  public:
    writer_state() { buffer_.reserve(buffer_reserve); }

    void append(const static_molecule& ligand) {
      buffer_ += ligand.properties.get(property_type::NAME);
      buffer_ += ' ';
      buffer_ += ligand.properties.get(property_type::SCORE);
      buffer_ += '\n';
      if (++lines_ % flush_every == 0) {
        flush();
      }
    }

    void flush() {
      std::cout << buffer_;
      buffer_.clear();
    }
  };

  class writer_body {
    writer_state& state_;
    global_ligands_limiter_t& global_ligands_limiter_;

  public:
    writer_body(writer_state& state, global_ligands_limiter_t& global_ligands_limiter)
        : state_(state), global_ligands_limiter_(global_ligands_limiter) {}

    tbb_flow::continue_msg operator()(lig_ptr ligand) const {
      state_.append(*ligand);
      global_ligands_limiter_.decrementer().try_put(tbb_flow::continue_msg{});
      return {};
    }
  };

} // namespace mudock::stream_engine
