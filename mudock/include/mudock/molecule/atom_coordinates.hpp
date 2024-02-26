#pragma once

#include <concepts>
#include <cstddef>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  template<class container_aliases>
  class atom_coordinates {
    template<typename T>
    using array_type = container_aliases::template atoms_size<T>;

    array_type<coordinate_type> x_coordinates;
    array_type<coordinate_type> y_coordinates;
    array_type<coordinate_type> z_coordinates;

  public:
    void resize(const std::size_t n);
    void fill(const coordinate_type value = 0);

    // utility functions to get the whole container
    [[nodiscard]] inline std::span<coordinate_type> x() { return x_coordinates; }
    [[nodiscard]] inline std::span<const coordinate_type> x() const { return x_coordinates; }
    [[nodiscard]] inline std::span<coordinate_type> y() { return y_coordinates; }
    [[nodiscard]] inline std::span<const coordinate_type> y() const { return y_coordinates; }
    [[nodiscard]] inline std::span<coordinate_type> z() { return z_coordinates; }
    [[nodiscard]] inline std::span<const coordinate_type> z() const { return z_coordinates; }

    // utility functions to access the data
    [[nodiscard]] inline coordinate_type& x(auto i) { return x_coordinates[i]; }
    [[nodiscard]] inline const coordinate_type& x(auto i) const { return x_coordinates[i]; }
    [[nodiscard]] inline coordinate_type& y(auto i) { return y_coordinates[i]; }
    [[nodiscard]] inline const coordinate_type& y(auto i) const { return y_coordinates[i]; }
    [[nodiscard]] inline coordinate_type& z(auto i) { return z_coordinates[i]; }
    [[nodiscard]] inline const coordinate_type& z(auto i) const { return z_coordinates[i]; }
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
