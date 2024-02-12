#pragma once

#include <mudock/type_alias.hpp>

namespace mudock {

  // this is the maximum number of atoms that we can expect from a static storage
  constexpr auto max_static_atoms() { return index_type{255}; }

} // namespace mudock
