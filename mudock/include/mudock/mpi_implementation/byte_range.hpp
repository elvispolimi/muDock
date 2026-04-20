#pragma once

#include <cstdint>

namespace mudock {

  struct byte_range {
    std::uint64_t begin = 0;
    std::uint64_t end   = 0;

    [[nodiscard]] bool empty() const { return begin >= end; }
  };

} // namespace mudock
