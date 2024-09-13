#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <mudock/molecule/constraints.hpp>
#include <span>
#include <vector>

namespace mudock {

  // this set of type alias describe a molecule that allocates a fixed size of memory to hold atoms and bonds
  // properties. On the bright side, we can reduce the number of allocs, but we need to introduce constraints
  // on the maximum number of atoms and bonds that we can process. We typically use this setup for ligands.
  struct static_containers {
    template<typename T>
    using atoms_size = std::array<T, max_static_atoms()>;

    template<typename T>
    using bonds_size = std::array<T, max_static_bonds()>;

    template<typename T>
    using fragments_size = std::array<T, max_static_atoms() * max_static_bonds()>;
  };

  // this set of type alas describe a molecule that allocates the required memory at runtime.It is the most
  // flexible approch, but it introduces runtime overheads. We typically use this setup for the proteins.
  struct dynamic_containers {
    template<typename T>
    using atoms_size = std::vector<T>;

    template<typename T>
    using bonds_size = std::vector<T>;

    template<typename T>
    using fragments_size = std::vector<T>;
  };

  // this concepts tell if a type is a specification of the container abstraction used in a class.
  template<class T>
  concept is_container_specification = (std::same_as<std::remove_cvref_t<T>, static_containers> ||
                                        std::same_as<std::remove_cvref_t<T>, dynamic_containers>);

  // helper functions to resize a container
  template<typename value_type>
  inline void resize(std::vector<value_type>& container, const int new_size) {
    container.resize(new_size);
  }
  template<typename value_type, std::size_t num_elements>
  inline void resize([[maybe_unused]] std::array<value_type, num_elements>& container,
                     [[maybe_unused]] const int new_size) {
    assert(static_cast<std::size_t>(new_size) <= num_elements);
  }

  // helper functions to fill a container with a value
  template<typename value_type>
  inline void fill(std::vector<value_type>& container, const value_type value) {
    std::fill(std::begin(container), std::end(container), value);
  }
  template<typename value_type, std::size_t num_elements>
  inline void fill(std::array<value_type, num_elements>& container, const value_type value) {
    container.fill(value);
  }

  // helper functions to remove an element by index
  template<typename value_type>
  inline void remove_atom(std::vector<value_type>& container, const int index) {
    assert(container.size() > index);
    std::shift_left(std::begin(container) + index, std::end(container), 1);
    container.resize(container.size() - std::size_t{1});
  }
  template<typename value_type, std::size_t num_elements>
  inline void remove_atom(std::array<value_type, num_elements>& container, const int index) {
    assert(container.size() > index);
    std::shift_left(std::begin(container) + index, std::end(container), 1);
  }

  // helper functions to create a span from a container
  template<typename value_type>
  [[nodiscard]] inline auto make_span(std::vector<value_type>& container, const int size) {
    return std::span(std::begin(container), size);
  }
  template<typename value_type, std::size_t num_elements>
  [[nodiscard]] inline auto make_span(std::array<value_type, num_elements>& container, const int size) {
    return std::span(std::begin(container), size);
  }
  template<typename value_type>
  [[nodiscard]] inline auto make_span(const std::vector<value_type>& container, const int size) {
    return std::span(std::begin(container), size);
  }
  template<typename value_type, std::size_t num_elements>
  [[nodiscard]] inline auto make_span(const std::array<value_type, num_elements>& container, const int size) {
    return std::span(std::begin(container), size);
  }

} // namespace mudock
