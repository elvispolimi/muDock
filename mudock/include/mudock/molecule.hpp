#pragma once

#include <cassert>
#include <concepts>
#include <cstdint>
#include <mudock/chem/autodock_types.hpp>
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
    requires is_container_specification<container_aliases>
  class molecule {
    template<typename T>
    using atoms_array_type = container_aliases::template atoms_size<T>;

    template<typename T>
    using bonds_array_type = container_aliases::template bonds_size<T>;

    // the atoms chemical properties
    atoms_array_type<element> atom_elements;
    // the atoms autodock types
    atoms_array_type<autodock_ff> atom_autodock_elements;
    std::size_t atom_size = std::size_t{0};

    // Atom "i" is aromatic, needed to get the autodock type
    atoms_array_type<int> is_aromatic;

    // the intra-molecular connections
    bonds_array_type<bond> bond_descriptions;
    std::size_t bonds_size = std::size_t{0};

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
    [[nodiscard]] inline auto elements() { return std::span(std::begin(atom_elements), atom_size); }
    [[nodiscard]] inline auto elements() const { return std::span(std::cbegin(atom_elements), atom_size); }
    [[nodiscard]] inline auto autodock_elements() {
      return std::span(std::begin(atom_autodock_elements), atom_size);
    }
    [[nodiscard]] inline auto autodock_elements() const {
      return std::span(std::cbegin(atom_autodock_elements), atom_size);
    }
    [[nodiscard]] inline auto aromaticity() { return std::span(std::begin(is_aromatic), atom_size); }
    [[nodiscard]] inline auto aromaticity() const { return std::span(std::cbegin(is_aromatic), atom_size); }
    [[nodiscard]] inline auto bonds() { return std::span(std::begin(bond_descriptions), bonds_size); }
    [[nodiscard]] inline auto bonds() const { return std::span(std::cbegin(bond_descriptions), bonds_size); }

    // TODO this invalidate other elements
    inline void remove_atom(const size_t index) {
      coordinates.remove_atom(index);
      atom_size--;

      atom_elements.erase(atom_elements.begin() + index);
      atom_autodock_elements.erase(atom_autodock_elements.begin() + index);
      is_aromatic.erase(is_aromatic.begin() + index);
      std::remove_if(bond_descriptions.begin(), bond_descriptions.end(), [&](auto bond) {
        return bond.source == index || bond.dest == index;
      });
      bonds_size = bond_descriptions.size();
    };

    // utility functions for accessing directly to a specific element
    [[nodiscard]] inline element& elements(std::size_t i) {
      assert(i < atom_size);
      return atom_elements[i];
    }
    [[nodiscard]] inline const element& elements(std::size_t i) const {
      assert(i < atom_size);
      return atom_elements[i];
    }
    [[nodiscard]] inline autodock_ff& autodock_elements(std::size_t i) {
      assert(i < atom_size);
      return atom_autodock_elements[i];
    }
    [[nodiscard]] inline const autodock_ff& autodock_elements(std::size_t i) const {
      assert(i < atom_size);
      return atom_autodock_elements[i];
    }
    [[nodiscard]] inline int& aromaticity(std::size_t i) {
      assert(i < atom_size);
      return is_aromatic[i];
    }
    [[nodiscard]] inline const int& aromaticity(std::size_t i) const {
      assert(i < atom_size);
      return is_aromatic[i];
    }
    [[nodiscard]] inline bond& bonds(std::size_t i) {
      assert(i < bonds_size);
      return bond_descriptions[i];
    }
    [[nodiscard]] inline const bond& bonds(std::size_t i) const {
      assert(i < bonds_size);
      return bond_descriptions[i];
    }
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
