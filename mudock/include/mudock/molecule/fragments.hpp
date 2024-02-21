#pragma once

#include <cstdint>
#include <mudock/grid/mdindex.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  template<template<typename> class container_type>
  class fragments {
    container_type<coordinate_type> storage;
    index2D index;

  public:
    void reset(const std::size_t num_atoms, const std::size_t num_bonds);

    inline std::span<coordinate_type> get_mask(const std::size_t bond_index) const {
      // NOTE: we assume that container_type will hold the elements linearly
      return {&storage[index.to1D(0, bond_index)], index.size_x()};
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
