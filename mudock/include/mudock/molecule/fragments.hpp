#pragma once

#include <concepts>
#include <cstdint>
#include <mudock/grid/mdindex.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<template<typename> class container_type>
  class fragments {
    container_type<coordinate_type> storage;
    index2D index;

  public:
    void reset(const std::size_t num_atoms, const std::size_t num_bonds);

    inline const coordinate_type* get_const_mask(const std::size_t bond_index) const {
      return &storage[index.to1D(0, bond_index)];
    }
    inline coordinate_type* get_mask(const std::size_t bond_index) const {
      return &storage[index.to1D(0, bond_index)];
    }
  };

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<>
  void fragments<static_container_type>::reset(const std::size_t num_atoms, const std::size_t num_bonds);

  template<>
  void fragments<dynamic_container_type>::reset(const std::size_t num_atoms, const std::size_t num_bonds);

} // namespace mudock
