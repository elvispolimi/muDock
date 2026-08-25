#pragma once

#include "stream_engine_types.hpp"

#include <atomic>
#include <cstddef>
#include <exception>
#include <memory>
#include <mudock/format/reader.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/molecule.hpp>
#include <oneapi/tbb/flow_graph.h>
#include <string>

namespace mudock::stream_engine {

  template<supported_format format>
  class parser_body {
    global_ligands_limiter_t& global_ligands_limiter_;
    std::atomic<std::size_t>* skipped_ligands_;
    std::atomic<bool>* stop_requested_;

  public:
    explicit parser_body(global_ligands_limiter_t& global_ligands_limiter,
                         std::atomic<std::size_t>* skipped = nullptr,
                         std::atomic<bool>* stop_signal    = nullptr)
        : global_ligands_limiter_(global_ligands_limiter),
          skipped_ligands_(skipped),
          stop_requested_(stop_signal) {}

    void operator()(const std::string& token, parser_mfn_t::output_ports_type& out) const {
      if (stop_requested_ != nullptr && stop_requested_->load(std::memory_order_relaxed)) {
        release_slot();
        return;
      }

      try {
        auto ligand = std::make_shared<static_molecule>(mudock::parser<format, static_molecule>(token));
        std::get<0>(out).try_put(std::move(ligand));
      } catch (const std::exception&) {
        if (skipped_ligands_ != nullptr) {
          skipped_ligands_->fetch_add(1, std::memory_order_relaxed);
        }
        release_slot();
      }
    }

  private:
    void release_slot() const { global_ligands_limiter_.decrementer().try_put(tbb_flow::continue_msg{}); }
  };

} // namespace mudock::stream_engine
