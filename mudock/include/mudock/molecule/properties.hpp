#pragma once

#include <cstddef>
#include <cstdint>

namespace mudock {
  enum class property_type : int { NAME = 0, SCORE, SOURCE_PATH, END };

  [[nodiscard]] constexpr int property_type_size() { return static_cast<int>(property_type::END); }
} // namespace mudock
