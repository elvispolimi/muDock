#pragma once

#include <mudock/compute/scratchpad_soa.hpp>

namespace mudock {

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  using scratchpad = scratchpad_soa<queue_type>;

}  // namespace mudock
