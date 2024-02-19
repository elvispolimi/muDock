#include <cassert>
#include <cstdint>
#include <mudock/grid.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/molecule/fragments.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<>
  void fragments<static_container_type>::reset(const std::size_t num_atoms, const std::size_t num_bonds) {
    assert(num_atoms < max_static_atoms());
    assert(num_bonds < max_static_bonds());
    storage.fill(coordinate_type{0});
    index = index2D{num_atoms, num_bonds};
  }

  template<>
  void fragments<dynamic_container_type>::reset(const std::size_t num_atoms, const std::size_t num_bonds) {
    storage.clear();
    storage.resize(num_atoms * num_bonds, coordinate_type{0});
    index = index2D{num_atoms, num_bonds};
  }

} // namespace mudock
