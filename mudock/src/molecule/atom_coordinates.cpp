#include <algorithm>
#include <cassert>
#include <mudock/molecule/atom_coordinates.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<>
  void atom_coordinates_type<static_container_type>::resize(const std::size_t n) {
    assert(n <= max_static_atoms());
  }

  template<>
  void atom_coordinates_type<dynamic_container_type>::resize(const std::size_t n) {
    x.resize(n, coordinate_type{0});
    y.resize(n, coordinate_type{0});
    z.resize(n, coordinate_type{0});
  }

  template<>
  void atom_coordinates_type<static_container_type>::fill(const coordinate_type value) {
    x.fill(coordinate_type{value});
    y.fill(coordinate_type{value});
    z.fill(coordinate_type{value});
  }

  template<>
  void atom_coordinates_type<dynamic_container_type>::fill(const coordinate_type value) {
    std::fill(std::begin(x), std::end(x), coordinate_type{value});
    std::fill(std::begin(y), std::end(y), coordinate_type{value});
    std::fill(std::begin(z), std::end(z), coordinate_type{value});
  }

} // namespace mudock
