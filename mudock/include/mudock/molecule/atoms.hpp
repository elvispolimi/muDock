#pragma once

#include <cstddef>
#include <mudock/chem/elements.hpp>
#include <mudock/molecule/atom_coordinates.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<template<typename> class container_type>
  class atoms_type {
  public:
    void resize(const std::size_t n);

    atom_coordinates_type<container_type> coordinates;
    container_type<element> elements;
    index_type num_atoms = index_type{0};
  };

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<>
  void atoms_type<static_container_type>::resize(const std::size_t n);
  template<>
  void atoms_type<dynamic_container_type>::resize(const std::size_t n);

} // namespace mudock
