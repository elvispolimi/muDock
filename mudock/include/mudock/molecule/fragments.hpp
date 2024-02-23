#pragma once

#include <cstdint>
#include <mudock/grid/mdindex.hpp>
#include <mudock/molecule/bond.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  template<template<typename> class container_type>
  class fragments {
    container_type<coordinate_type> storage;
    index2D index;

  public:
    fragments(const std::size_t num_atoms, const std::size_t num_bonds);

    inline std::span<coordinate_type> get_mask(const std::size_t bond_index) {
      const auto begin = storage.begin() + bond_index;
      const auto end   = begin + index.size_x();
      assert(begin != std::end(storage) && end != std::end(storage));
      return std::span(begin, end);
    }
    inline std::span<const coordinate_type> get_mask(const std::size_t bond_index) const {
      const auto begin = storage.cbegin() + bond_index;
      const auto end   = begin + index.size_x();
      assert(begin != std::end(storage) && end != std::end(storage));
      return std::span(begin, end);
    }
    inline auto get_num_rotatable_bonds() const { return index.size_y(); }
  };

  template<template<typename> class container_type>
  fragments<container_type> make_fragments(const container_type<bond>& bonds,
                                           const index_type num_atoms,
                                           const index_type num_bonds);

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<>
  fragments<static_container_type>::fragments(const std::size_t num_atoms, const std::size_t num_bonds);

  template<>
  fragments<dynamic_container_type>::fragments(const std::size_t num_atoms, const std::size_t num_bonds);

  template<>
  fragments<static_container_type>
      make_fragments<static_container_type>(const static_container_type<bond>& bonds,
                                            const index_type num_atoms,
                                            const index_type num_bonds);
  template<>
  fragments<dynamic_container_type>
      make_fragments<dynamic_container_type>(const dynamic_container_type<bond>& bonds,
                                             const index_type num_atoms,
                                             const index_type num_bonds);

} // namespace mudock
