#pragma once

#include <concepts>
#include <mudock/molecule/atoms.hpp>
#include <mudock/molecule/containers.hpp>
#include <string>

namespace mudock {

  // this is the generic definition of a molecule, that depends on the used type of storage
  template<template<typename> class container_type>
  struct molecule_type {
    std::string name;
    atoms_type<container_type> atoms;
  };

  // this is the concept that defines a molecule, which is any molecule_type for which we have defined a
  // special container and we are agnostic about it.
  template<class T>
  concept molecule = (std::same_as<T, molecule_type<static_container_type>> ||
                      std::same_as<T, molecule_type<dynamic_container_type>>);

} // namespace mudock
