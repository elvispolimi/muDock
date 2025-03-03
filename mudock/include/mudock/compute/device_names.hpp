#pragma once

#include <string_view>

namespace mudock {
  static constexpr auto cpu_token = std::string_view{"CPU"};
  static constexpr auto gpu_token = std::string_view{"GPU"};
} // namespace mudock
