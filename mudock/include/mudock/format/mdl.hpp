#pragma once

#include <string_view>

namespace mudock {
  struct sdf {
    std::string_view::size_type next_molecule_start_index(std::string_view text) const {
      constexpr std::string_view delimiter = "$$$$";
      const auto start = text.find(delimiter);
      if (start == std::string_view::npos) return start;
      const auto newline = text.find('\n', start + delimiter.size());
      if (newline == std::string_view::npos ||
          text.find_first_not_of(" \t\r\n", newline + 1) == std::string_view::npos)
        return std::string_view::npos;
      return newline + 1;
    }
  };

  struct mol {
    std::string_view::size_type next_molecule_start_index(std::string_view) const {
      return std::string_view::npos;
    }
  };
} // namespace mudock
