#include <algorithm>
#include <cassert>
#include <mudock/molecule/atom_coordinates.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<>
  void atom_coordinates<static_containers>::resize([[maybe_unused]] const std::size_t n) {
    assert(n <= max_static_atoms());
  }

  template<>
  void atom_coordinates<dynamic_containers>::resize(const std::size_t n) {
    x_coordinates.resize(n, coordinate_type{0});
    y_coordinates.resize(n, coordinate_type{0});
    z_coordinates.resize(n, coordinate_type{0});
  }

  template<>
  void atom_coordinates<static_containers>::fill(const coordinate_type value) {
    x_coordinates.fill(coordinate_type{value});
    y_coordinates.fill(coordinate_type{value});
    z_coordinates.fill(coordinate_type{value});
  }

  template<>
  void atom_coordinates<dynamic_containers>::fill(const coordinate_type value) {
    std::fill(std::begin(x_coordinates), std::end(x_coordinates), coordinate_type{value});
    std::fill(std::begin(y_coordinates), std::end(y_coordinates), coordinate_type{value});
    std::fill(std::begin(z_coordinates), std::end(z_coordinates), coordinate_type{value});
  }

} // namespace mudock
