#pragma once

#include <cstdint>
#include <gsl/pointers>
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

    inline gsl::not_null<coordinate_type*> get_mask(const std::size_t bond_index) {
      // NOTE: we assume that container_type will hold the elements linearly
      return &storage[index.to1D(0, bond_index)];
    }
    inline gsl::not_null<const coordinate_type*> get_mask(const std::size_t bond_index) const {
      // NOTE: we assume that container_type will hold the elements linearly
      return &storage[index.to1D(0, bond_index)];
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
