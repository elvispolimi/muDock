#pragma once

#include <cstddef>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<template<typename> class container_type>
  class atom_coordinates_type {
  public:
    void resize(const std::size_t n);
    void fill(const coordinate_type value = 0);

    container_type<coordinate_type> x;
    container_type<coordinate_type> y;
    container_type<coordinate_type> z;
  };

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<>
  void atom_coordinates_type<static_container_type>::resize(const std::size_t n);
  template<>
  void atom_coordinates_type<dynamic_container_type>::resize(const std::size_t n);

  template<>
  void atom_coordinates_type<static_container_type>::fill(const coordinate_type value);
  template<>
  void atom_coordinates_type<dynamic_container_type>::fill(const coordinate_type value);

} // namespace mudock
