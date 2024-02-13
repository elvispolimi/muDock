#include "mudock/chem/elements.hpp"
#include "mudock/type_alias.hpp"

#include <cassert>
#include <mudock/molecule/atoms.hpp>

namespace mudock {

  template<>
  void atoms_type<static_container_type>::resize(const std::size_t n) {
    coordinates.resize(n);
    assert(n <= max_static_atoms());
    num_atoms = static_cast<index_type>(n);
  }

  template<>
  void atoms_type<dynamic_container_type>::resize(const std::size_t n) {
    coordinates.resize(n);
    elements.resize(n, element::H);
    num_atoms = static_cast<index_type>(n);
  }

} // namespace mudock
