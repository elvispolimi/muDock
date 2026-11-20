#pragma once

#include "stream_engine_types.hpp"

#include <functional>
#include <memory>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/molecule.hpp>
#include <mudock/tbb_implementation/run_metrics.hpp>
#include <optional>
#include <utility>

namespace mudock::stream_engine {

  /**
   * Owns the reorder-buffer and shared_ptr-wraps batches as they fill.
   */
  class batcher_state {
    reorder_buffer<static_molecule> rob_;

  public:
    explicit batcher_state(std::function<int(const int)> get_size): rob_(std::move(get_size)) {}

    std::optional<lig_batch_ptr> add(lig_ptr ligand) {
      auto uptr             = std::make_unique<static_molecule>(std::move(*ligand));
      auto [batch, is_full] = rob_.add_ligand(std::move(uptr));
      if (is_full) {
        return std::optional{std::make_shared<lig_batch>(std::move(batch))};
      }
      return std::nullopt;
    }

    template<typename Sink>
    void flush_all(Sink sink) {
      while (true) {
        auto [batch, has_batches] = rob_.flush_one();
        if (!has_batches) {
          break;
        }
        sink(std::make_shared<lig_batch>(std::move(batch)));
      }
    }
  };

  /**
   * Multifunction-node body that forwards each ligand to the `batcher_state`
   * and emits a batch on port 0 when one fills.
   */
  class batcher_body {
    batcher_state& state_;
    observability::backend_metrics& metrics_;

  public:
    explicit batcher_body(batcher_state& state, observability::backend_metrics& metrics)
        : state_(state), metrics_(metrics) {}

    void operator()(lig_ptr ligand, batcher_mfn::output_ports_type& out) const {
      metrics_.ligands_batched.increase();
      if (const auto batch = state_.add(std::move(ligand))) {
        metrics_.batches_emitted.increase();
        std::get<0>(out).try_put(*batch);
      }
    }
  };

} // namespace mudock::stream_engine
