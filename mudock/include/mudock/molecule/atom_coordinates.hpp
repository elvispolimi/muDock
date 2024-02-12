#pragma once

#include <cstddef>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  /**
   * This struct represents all the atoms coordinates. All of them are zero-initialized. We expect that
   * external functions will fill with the actual content.
  */
  template<template<typename> class container_type>
  struct atom_coordinates_type {
    atom_coordinates_type(const std::size_t n);
    atom_coordinates_type()                                        = default;
    atom_coordinates_type(atom_coordinates_type&&)                 = default;
    atom_coordinates_type(const atom_coordinates_type&)            = default;
    ~atom_coordinates_type()                                       = default;
    atom_coordinates_type& operator=(atom_coordinates_type&&)      = default;
    atom_coordinates_type& operator=(const atom_coordinates_type&) = default;

    container_type<coordinate_type> x;
    container_type<coordinate_type> y;
    container_type<coordinate_type> z;
  };

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<>
  atom_coordinates_type<static_container_type>::atom_coordinates_type(const std::size_t n);
  template<>
  atom_coordinates_type<dynamic_container_type>::atom_coordinates_type(const std::size_t n);

} // namespace mudock
