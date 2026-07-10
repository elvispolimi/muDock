#pragma once

#include "stream_engine_types.hpp"

#include <memory>
#include <mudock/molecule.hpp>
#include <mudock/molecule/properties.hpp>
#include <utility>

namespace mudock::stream_engine {

  class scorer_body {
  public:
    void operator()(lig_batch_ptr batch, scorer_mfn_t::output_ports_type& out) const {
      for (int i = 0; i < batch->num_ligands; ++i) {
        auto& lig = batch->molecules[i];
        if (!lig) {
          continue;
        }
        lig->properties.assign(property_type::SCORE, "0.0");
        std::get<0>(out).try_put(lig_ptr{std::move(lig)});
      }
    }
  };

} // namespace mudock::stream_engine
