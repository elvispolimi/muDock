#pragma once

#include <cassert>
#include <concepts>
#include <cstddef>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  template<class container_aliases>
    requires is_container_specification<container_aliases>
  class atom_coordinates {
    template<typename T>
    using array_type = container_aliases::template atoms_size<T>;

    std::size_t num_atoms;
    array_type<coordinate_type> x_coordinates;
    array_type<coordinate_type> y_coordinates;
    array_type<coordinate_type> z_coordinates;

  public:
    void resize(const std::size_t n);
    void fill(const coordinate_type value = 0);

    // utility functions to get the whole container
    [[nodiscard]] inline auto x() { return std::span(std::begin(x_coordinates), num_atoms); }
    [[nodiscard]] inline auto x() const { return std::span(std::cbegin(x_coordinates), num_atoms); }
    [[nodiscard]] inline auto y() { return std::span(std::begin(y_coordinates), num_atoms); }
    [[nodiscard]] inline auto y() const { return std::span(std::cbegin(y_coordinates), num_atoms); }
    [[nodiscard]] inline auto z() { return std::span(std::begin(z_coordinates), num_atoms); }
    [[nodiscard]] inline auto z() const { return std::span(std::cbegin(z_coordinates), num_atoms); }

    // utility functions to access the data
    [[nodiscard]] inline coordinate_type& x(const std::size_t i) {
      assert(i < num_atoms);
      return x_coordinates[i];
    }
    [[nodiscard]] inline const coordinate_type& x(const std::size_t i) const {
      assert(i < num_atoms);
      return x_coordinates[i];
    }
    [[nodiscard]] inline coordinate_type& y(const std::size_t i) {
      assert(i < num_atoms);
      return y_coordinates[i];
    }
    [[nodiscard]] inline const coordinate_type& y(const std::size_t i) const {
      assert(i < num_atoms);
      return y_coordinates[i];
    }
    [[nodiscard]] inline coordinate_type& z(const std::size_t i) {
      assert(i < num_atoms);
      return z_coordinates[i];
    }
    [[nodiscard]] inline const coordinate_type& z(const std::size_t i) const {
      assert(i < num_atoms);
      return z_coordinates[i];
    }
  };

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<>
  void atom_coordinates<static_containers>::resize(const std::size_t n);
  template<>
  void atom_coordinates<dynamic_containers>::resize(const std::size_t n);

  template<>
  void atom_coordinates<static_containers>::fill(const coordinate_type value);
  template<>
  void atom_coordinates<dynamic_containers>::fill(const coordinate_type value);

} // namespace mudock
