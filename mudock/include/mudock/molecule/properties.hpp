#pragma once

#include <cstddef>
#include <cstdint>

namespace mudock {
  enum class property_type : std::size_t { NAME = 0, SCORE, END };

  [[nodiscard]] constexpr std::size_t property_type_size() {
    return static_cast<std::size_t>(property_type::END);
  }
} // namespace mudock
