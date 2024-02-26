#pragma once

#include <concepts>
#include <mudock/chem/elements.hpp>
#include <mudock/molecule/atom_coordinates.hpp>
#include <mudock/molecule/bond.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/molecule/fragments.hpp>
#include <mudock/molecule/graph.hpp>
#include <mudock/molecule/properties.hpp>
#include <mudock/molecule/property_table.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  // this is the generic definition of a molecule, that depends on the used type of storage
  template<class container_aliases>
  class molecule {
    template<typename T>
    using atoms_array_type = container_aliases::template atoms_size<T>;

    template<typename T>
    using bonds_array_type = container_aliases::template bonds_size<T>;

    // the atoms chemical properties
    atoms_array_type<element> atom_elements;
    index_type atom_size = index_type{0};

    // the intra-molecular connections
    bonds_array_type<bond> bond_descriptions;
    index_type bonds_size = index_type{0};

  public:
    // utility function to allocate the memory based on the number of atoms and bonds
    void resize(const std::size_t n_atoms, std::size_t n_bonds);

    // a container that we can use to store key-value properties, e.g. its name
    property_map properties;

    // the atoms coordinates
    atom_coordinates<container_aliases> coordinates;

    // utility functions to get the molecule geometry
    [[nodiscard]] constexpr auto num_atoms() const { return atom_size; }
    [[nodiscard]] constexpr auto num_bonds() const { return bonds_size; }

    // utility functions to get the whole containers for each fields
    [[nodiscard]] std::span<element> elements() { return atom_elements; }
    [[nodiscard]] std::span<const element> elements() const { return atom_elements; }
    [[nodiscard]] std::span<bond> bonds() { return bond_descriptions; }
    [[nodiscard]] std::span<const bond> bonds() const { return bond_descriptions; }

    // utility functions for accessing directly to a specific element
    [[nodiscard]] inline element& elements(auto i) { return atom_elements[i]; }
    [[nodiscard]] inline const element& elements(auto i) const { return atom_elements[i]; }
    [[nodiscard]] inline bond& bonds(auto i) { return bond_descriptions[i]; }
    [[nodiscard]] inline const bond& bonds(auto i) const { return bond_descriptions[i]; }
  };

  //===------------------------------------------------------------------------------------------------------
  // Type alias to deal with concrete molecule type (according to the storage type)
  //===------------------------------------------------------------------------------------------------------

  using dynamic_molecule = molecule<dynamic_containers>;
  using static_molecule  = molecule<static_containers>;

  // this is the concept that defines a molecule, which is any molecule for which we have defined a
  // special container and we are agnostic about it.
  template<class T>
  concept is_molecule = (std::same_as<std::remove_cvref_t<T>, static_molecule> ||
                         std::same_as<std::remove_cvref_t<T>, dynamic_molecule>);

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<>
  void molecule<dynamic_containers>::resize(const std::size_t n_atoms, std::size_t n_bonds);
  template<>
  void molecule<static_containers>::resize(const std::size_t n_atoms, std::size_t n_bonds);

} // namespace mudock
