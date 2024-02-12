#include <cassert>
#include <mudock/molecule/atom_coordinates.hpp>

namespace mudock {

  template<>
  atom_coordinates_type<static_container_type>::atom_coordinates_type(const std::size_t n) {
    assert(n <= max_static_atoms());
    x.fill(coordinate_type{0});
    y.fill(coordinate_type{0});
    z.fill(coordinate_type{0});
  }

  template<>
  atom_coordinates_type<dynamic_container_type>::atom_coordinates_type(const std::size_t n)
      : x(n, coordinate_type{0}), y(n, coordinate_type{0}), z(n, coordinate_type{0}) {}

} // namespace mudock
