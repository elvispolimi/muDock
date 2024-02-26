#pragma once

#include <cstdint>

namespace mudock {

  // this is the maximum number of atoms that we can expect from a static storage
  [[nodiscard]] constexpr auto max_static_atoms() { return std::size_t{255}; }

  // this is the maximum number of bonds that we can expect from a static storage
  [[nodiscard]] constexpr auto max_static_bonds() { return max_static_atoms() * std::size_t{2}; }

} // namespace mudock
