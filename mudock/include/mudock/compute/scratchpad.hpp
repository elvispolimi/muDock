#pragma once

#include <mudock/compute/scratchpad_soa.hpp>
#include <mudock/compute/scratchpad_layout_sets.hpp>

namespace mudock {

  template<typename queue_type, typename layout_set = all_scratch_layouts>
    requires std::derived_from<queue_type, queue>
  using scratchpad = scratchpad_soa<queue_type, layout_set>;

}  // namespace mudock
