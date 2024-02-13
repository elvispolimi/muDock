#pragma once

#include <concepts>
#include <mudock/molecule/atoms.hpp>
#include <mudock/molecule/bond.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <string>

namespace mudock {

  // this is the generic definition of a molecule, that depends on the used type of storage
  template<template<typename> class container_type>
  class molecule_type {
  public:
    void resize(const std::size_t n_atoms, std::size_t n_bonds);

    std::string name;
    atoms_type<container_type> atoms;
    container_type<bond> bonds;
    index_type num_bonds = index_type{0};
  };
  template<>
  void molecule_type<static_container_type>::resize(const std::size_t n_atoms, std::size_t n_bonds);
  template<>
  void molecule_type<dynamic_container_type>::resize(const std::size_t n_atoms, std::size_t n_bonds);

  using dynamic_molecule = molecule_type<dynamic_container_type>;
  using static_molecule  = molecule_type<static_container_type>;

  // this is the concept that defines a molecule, which is any molecule_type for which we have defined a
  // special container and we are agnostic about it.
  template<class T>
  concept molecule = (std::same_as<T, static_molecule> || std::same_as<T, dynamic_molecule>);

} // namespace mudock
