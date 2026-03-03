#pragma once

#include <algorithm>

namespace mudock {
  struct batch_multiple {
    int active_blocks_per_sm{1};
    int num_sms{1};

    int total_multiple() const {
      return std::max(1, active_blocks_per_sm) * std::max(1, num_sms);
    }
  };

  inline batch_multiple normalize_batch_multiple(batch_multiple value) {
    value.active_blocks_per_sm = std::max(1, value.active_blocks_per_sm);
    value.num_sms              = std::max(1, value.num_sms);
    return value;
  }
} // namespace mudock
