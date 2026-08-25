#pragma once

#include "stream_engine_types.hpp"

#include <memory>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/molecule.hpp>
#include <optional>
#include <utility>

namespace mudock::stream_engine {

  class batcher_state {
    reorder_buffer<static_molecule> rob_;

  public:
    std::optional<lig_batch_t> add(const lig_ptr& ligand) {
      auto uptr             = std::make_unique<static_molecule>(std::move(*ligand));
      auto [batch, is_full] = rob_.add_ligand(std::move(uptr));
      if (is_full) {
        return std::optional{std::move(batch)};
      }
      return std::nullopt;
    }

    template<typename Sink>
    void flush_all(Sink&& sink) {
      while (true) {
        auto [batch, has_batches] = rob_.flush_one();
        if (!has_batches) {
          break;
        }
        sink(std::move(batch));
      }
    }
  };

  class batcher_body {
    batcher_state& state_;

  public:
    explicit batcher_body(batcher_state& state): state_(state) {}

    void operator()(lig_ptr ligand, batcher_mfn_t::output_ports_type& out) const {
      if (auto batch = state_.add(std::move(ligand))) {
        std::get<0>(out).try_put(std::make_shared<lig_batch_t>(std::move(*batch)));
      }
    }
  };

} // namespace mudock::stream_engine
