#pragma once

#include <mudock/type_alias.hpp>

namespace mudock {

  // this is the maximum number of atoms that we can expect from a static storage
  constexpr auto max_static_atoms() { return index_type{255}; }

  // this is the maximum number of bonds that we can expect from a static storage
  constexpr auto max_static_bonds() { return max_static_atoms() * index_type{2}; }

} // namespace mudock
