#pragma once

#include <cstddef>
#include <mudock/molecule/atom_coordinates.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  /**
   * This struct represents all the properties that we use to describe the molecule atoms. All the properties
   * are zero-initialized. We expect that external functions will fill with the actual content.
  */
  template<template<typename> class container_type>
  struct atoms_type {
    atoms_type(): num_atoms(index_type{0}) {}
    atoms_type(const std::size_t n);
    atoms_type(atoms_type&&)                 = default;
    atoms_type(const atoms_type&)            = default;
    ~atoms_type()                            = default;
    atoms_type& operator=(atoms_type&&)      = default;
    atoms_type& operator=(const atoms_type&) = default;

    atom_coordinates_type<container_type> coordinates;
    index_type num_atoms;
  };

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<>
  atoms_type<static_container_type>::atoms_type(const std::size_t n);
  template<>
  atoms_type<dynamic_container_type>::atoms_type(const std::size_t n);

} // namespace mudock
