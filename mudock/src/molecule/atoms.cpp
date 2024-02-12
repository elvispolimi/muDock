#include <cassert>
#include <mudock/molecule/atoms.hpp>

namespace mudock {

  template<>
  atoms_type<static_container_type>::atoms_type(const std::size_t n)
      : coordinates(n), num_atoms(static_cast<index_type>(n)) {}

  template<>
  atoms_type<dynamic_container_type>::atoms_type(const std::size_t n)
      : coordinates(n), num_atoms(static_cast<index_type>(n)) {}

} // namespace mudock
