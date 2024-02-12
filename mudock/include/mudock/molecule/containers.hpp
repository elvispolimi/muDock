#pragma once

#include <array>
#include <mudock/molecule/constraints.hpp>
#include <vector>

namespace mudock {

  // this container does not use runtime memory, but it is equivalent to a plain old array. For this reason
  // the memory will be allocated depending on the variable definition. The drawback that we need to tell
  // a maximum number of atoms that we can support.
  template<typename value_type>
  using static_container_type = std::array<value_type, max_static_atoms()>;

  // this container allocate memory at runtime depending on the actual requirements. It's a generalization of
  // the previous one, but it has higher overhead.
  template<typename value_type>
  using dynamic_container_type = std::vector<value_type>;

} // namespace mudock
