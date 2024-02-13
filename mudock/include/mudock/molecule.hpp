#pragma once

#include <concepts>
#include <mudock/chem/elements.hpp>
#include <mudock/molecule/atom_coordinates.hpp>
#include <mudock/molecule/bond.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/molecule/property_table.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  // this is the generic definition of a molecule, that depends on the used type of storage
  template<template<typename> class container_type>
  class molecule {
  public:
    // the atoms description
    atom_coordinates<container_type> coordinates;
    container_type<element> elements;
    index_type num_atoms = index_type{0};

    // the bonds description
    container_type<bond> bonds;
    index_type num_bonds = index_type{0};

    // a container that we can use to store key-value properties, e.g. its name
    property_map properties;

    // utility function to allocate the memory based on the number of atoms and bonds
    void resize(const std::size_t n_atoms, std::size_t n_bonds);
  };

  //===------------------------------------------------------------------------------------------------------
  // Type alias to deal with concrete molecule type (according to the storage type)
  //===------------------------------------------------------------------------------------------------------

  using dynamic_molecule = molecule<dynamic_container_type>;
  using static_molecule  = molecule<static_container_type>;

  // this is the concept that defines a molecule, which is any molecule for which we have defined a
  // special container and we are agnostic about it.
  template<class T>
  concept is_molecule = (std::same_as<T, static_molecule> || std::same_as<T, dynamic_molecule>);

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<>
  void molecule<static_container_type>::resize(const std::size_t n_atoms, std::size_t n_bonds);
  template<>
  void molecule<dynamic_container_type>::resize(const std::size_t n_atoms, std::size_t n_bonds);

} // namespace mudock
