#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <mudock/controllers.hpp>
#include <mudock/controllers/squeue.hpp>
#include <mudock/devices.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  template<language_types l_t>
  void worker(std::shared_ptr<squeue<static_molecule>> i_queue,
              std::shared_ptr<squeue<static_molecule>> o_queue,
              const index_type id) {
    molecules_scratchpad<l_t> m_scratch;
    device_context<l_t> context{id};

    auto molecules = i_queue->dequeue(get_max_work());
    while (!molecules.empty()) {
      m_scratch.copy_in(molecules, context);
      dock_and_score(m_scratch, context);
      m_scratch.copy_out(molecules, context);
      o_queue->enqueue(molecules);
      molecules = i_queue->dequeue(get_max_work());
    }
  };
} // namespace mudock
