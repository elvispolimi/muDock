#pragma once

#include "stream_engine_types.hpp"

#include <utility>

namespace mudock::stream_engine {

  /**
   * Non-template base so a `backend_subgraph` can own a heterogeneous set of `screener_state<stage_t>` instances
   * behind a single vector without knowing `stage_t`.
   */
  struct screener_state_base {
    virtual ~screener_state_base() = default;
  };

  /**
   * Owns one stage instance (not thread-safe) and runs one batch through prepare / compute / teardown.
   */
  template<typename stage_t>
  class screener_state: public screener_state_base {
    stage_t stage_;

  public:
    explicit screener_state(stage_t stage): stage_(std::move(stage)) {}

    void screen(lig_batch& batch) {
      stage_.prepare(batch);
      stage_();
      stage_.teardown(batch);
    }
  };

  /**
   * Multifunction-node body that screens one batch via its `screener_state`
   * and emits each resulting ligand on port 0.
   */
  template<typename stage_t>
  class screener_body {
    screener_state<stage_t>& state_;

  public:
    explicit screener_body(screener_state<stage_t>& state): state_(state) {}

    void operator()(lig_batch_ptr batch, screener_mfn::output_ports_type& out) const {
      state_.screen(*batch);
      for (size_t i = 0; i < static_cast<size_t>(batch->num_ligands); ++i) {
        auto& lig = batch->molecules[i];
        if (!lig) {
          continue;
        }
        std::get<0>(out).try_put(lig_ptr{std::move(lig)});
      }
    }
  };

} // namespace mudock::stream_engine
