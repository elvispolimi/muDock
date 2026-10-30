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

  // batch_multiple.hpp
  inline batch_multiple min_batch_multiple(const batch_multiple a, const batch_multiple b) {
    const auto na = normalize_batch_multiple(a);
    const auto nb = normalize_batch_multiple(b);
    return (na.total_multiple() <= nb.total_multiple()) ? na : nb;
  }
} // namespace mudock
