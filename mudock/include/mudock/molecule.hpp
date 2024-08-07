#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <mudock/chem/elements.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/molecule/bond.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/molecule/containers.hpp>
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
  public:
    template<typename T>
    using atoms_array_type = container_aliases::template atoms_size<T>;

    template<typename T>
    using bonds_array_type = container_aliases::template bonds_size<T>;

  private:
    // the atoms chemical properties
    atoms_array_type<element> atom_elements;
    atoms_array_type<autodock_ff> atom_autodock_type;
    atoms_array_type<fp_type> x_coordinates;
    atoms_array_type<fp_type> y_coordinates;
    atoms_array_type<fp_type> z_coordinates;
    atoms_array_type<int> atom_is_aromatic; // 1-> aromatic, 0 -> no aromatic
    atoms_array_type<fp_type> atom_Rii;
    atoms_array_type<fp_type> atom_vol;
    atoms_array_type<fp_type> atom_solpar;
    atoms_array_type<fp_type> atom_epsii;
    atoms_array_type<fp_type> atom_Rij_hb;
    atoms_array_type<fp_type> atom_epsij_hb;
    atoms_array_type<fp_type> atom_charge;
    atoms_array_type<std::size_t> atom_num_hbond;
    std::size_t atoms_size = std::size_t{0};

    // the intra-molecular connections
    bonds_array_type<bond> bond_descriptions;
    std::size_t bonds_size = std::size_t{0};

  public:
    // functions to manage the geometry of a molecule
    void resize(const std::size_t n_atoms, std::size_t n_bonds);
    void remove_atom(const std::size_t index);
    [[nodiscard]] constexpr auto num_atoms() const { return atoms_size; }
    [[nodiscard]] constexpr auto num_bonds() const { return bonds_size; }
    constexpr auto num_rotamers() const {
      return std::count_if(std::begin(bond_descriptions), std::end(bond_descriptions), [](const bond& b) {
        return b.can_rotate;
      });
    }

    // a container that we can use to store key-value properties, e.g. its name
    property_map properties;

    // utility functions to get information about the bonds
    [[nodiscard]] inline auto get_bonds() { return std::span(std::begin(bond_descriptions), bonds_size); }
    [[nodiscard]] inline auto get_bonds() const {
      return std::span(std::cbegin(bond_descriptions), bonds_size);
    }
    [[nodiscard]] inline auto& bonds(std::size_t index) { return bond_descriptions[index]; }
    [[nodiscard]] inline const auto& bonds(std::size_t index) const { return bond_descriptions[index]; }

    // utility functions to get the span of the whole molecule (read + write)
    [[nodiscard]] inline auto get_elements() { return make_span(atom_elements, atoms_size); }
    [[nodiscard]] inline auto get_autodock_type() { return make_span(atom_autodock_type, atoms_size); }
    [[nodiscard]] inline auto get_x() { return make_span(x_coordinates, atoms_size); }
    [[nodiscard]] inline auto get_y() { return make_span(y_coordinates, atoms_size); }
    [[nodiscard]] inline auto get_z() { return make_span(z_coordinates, atoms_size); }
    [[nodiscard]] inline auto get_is_aromatic() { return make_span(atom_is_aromatic, atoms_size); }
    [[nodiscard]] inline auto get_Rii() { return make_span(atom_Rii, atoms_size); }
    [[nodiscard]] inline auto get_vol() { return make_span(atom_vol, atoms_size); }
    [[nodiscard]] inline auto get_solpar() { return make_span(atom_solpar, atoms_size); }
    [[nodiscard]] inline auto get_epsii() { return make_span(atom_epsii, atoms_size); }
    [[nodiscard]] inline auto get_Rij_hb() { return make_span(atom_Rij_hb, atoms_size); }
    [[nodiscard]] inline auto get_epsij_hb() { return make_span(atom_epsij_hb, atoms_size); }
    [[nodiscard]] inline auto get_charge() { return make_span(atom_charge, atoms_size); }
    [[nodiscard]] inline auto get_num_hbond() { return make_span(atom_num_hbond, atoms_size); }

    // utility functions to get the span of the whole molecule (read only)
    [[nodiscard]] inline auto get_element() const { return make_span(atom_elements, atoms_size); }
    [[nodiscard]] inline auto get_autodock_type() const { return make_span(atom_autodock_type, atoms_size); }
    [[nodiscard]] inline auto get_x() const { return make_span(x_coordinates, atoms_size); }
    [[nodiscard]] inline auto get_y() const { return make_span(y_coordinates, atoms_size); }
    [[nodiscard]] inline auto get_z() const { return make_span(z_coordinates, atoms_size); }
    [[nodiscard]] inline auto get_is_aromatic() const { return make_span(atom_is_aromatic, atoms_size); }
    [[nodiscard]] inline auto get_Rii() const { return make_span(atom_Rii, atoms_size); }
    [[nodiscard]] inline auto get_vol() const { return make_span(atom_vol, atoms_size); }
    [[nodiscard]] inline auto get_solpar() const { return make_span(atom_solpar, atoms_size); }
    [[nodiscard]] inline auto get_epsii() const { return make_span(atom_epsii, atoms_size); }
    [[nodiscard]] inline auto get_Rij_hb() const { return make_span(atom_Rij_hb, atoms_size); }
    [[nodiscard]] inline auto get_epsij_hb() const { return make_span(atom_epsij_hb, atoms_size); }
    [[nodiscard]] inline auto get_charge() const { return make_span(atom_charge, atoms_size); }
    [[nodiscard]] inline auto get_num_hbond() const { return make_span(atom_num_hbond, atoms_size); }

    // utility functions to get the ref to an atom element (read + write)
    [[nodiscard]] inline auto& elements(const std::size_t index) { return atom_elements[index]; }
    [[nodiscard]] inline auto& autodock_type(const std::size_t index) { return atom_autodock_type[index]; }
    [[nodiscard]] inline auto& x(const std::size_t index) { return x_coordinates[index]; }
    [[nodiscard]] inline auto& y(const std::size_t index) { return y_coordinates[index]; }
    [[nodiscard]] inline auto& z(const std::size_t index) { return z_coordinates[index]; }
    [[nodiscard]] inline auto& is_aromatic(const std::size_t index) { return atom_is_aromatic[index]; }
    [[nodiscard]] inline auto& Rii(const std::size_t index) { return atom_Rii[index]; }
    [[nodiscard]] inline auto& vol(const std::size_t index) { return atom_vol[index]; }
    [[nodiscard]] inline auto& solpar(const std::size_t index) { return atom_solpar[index]; }
    [[nodiscard]] inline auto& epsii(const std::size_t index) { return atom_epsii[index]; }
    [[nodiscard]] inline auto& Rij_hb(const std::size_t index) { return atom_Rij_hb[index]; }
    [[nodiscard]] inline auto& epsij_hb(const std::size_t index) { return atom_epsij_hb[index]; }
    [[nodiscard]] inline auto& charge(const std::size_t index) { return atom_charge[index]; }
    [[nodiscard]] inline auto& num_hbond(const std::size_t index) { return atom_num_hbond[index]; }

    // utility functions to get the span of the whole molecule (read only)
    [[nodiscard]] inline const auto& elements(const std::size_t index) const { return atom_elements[index]; }
    [[nodiscard]] inline const auto& autodock_type(const std::size_t index) const { return atom_autodock_type[index]; }
    [[nodiscard]] inline const auto& x(const std::size_t index) const { return x_coordinates[index]; }
    [[nodiscard]] inline const auto& y(const std::size_t index) const { return y_coordinates[index]; }
    [[nodiscard]] inline const auto& z(const std::size_t index) const { return z_coordinates[index]; }
    [[nodiscard]] inline const auto& is_aromatic(const std::size_t index) const {
      return atom_is_aromatic[index];
    }
    [[nodiscard]] inline const auto& Rii(const std::size_t index) const { return atom_Rii[index]; }
    [[nodiscard]] inline const auto& vol(const std::size_t index) const { return atom_vol[index]; }
    [[nodiscard]] inline const auto& solpar(const std::size_t index) const { return atom_solpar[index]; }
    [[nodiscard]] inline const auto& epsii(const std::size_t index) const { return atom_epsii[index]; }
    [[nodiscard]] inline const auto& Rij_hb(const std::size_t index) const { return atom_Rij_hb[index]; }
    [[nodiscard]] inline const auto& epsij_hb(const std::size_t index) const { return atom_epsij_hb[index]; }
    [[nodiscard]] inline auto& charge(const std::size_t index) const { return atom_charge[index]; }
    [[nodiscard]] inline const auto& num_hbond(const std::size_t index) const {
      return atom_num_hbond[index];
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

  template<class container_aliases>
    requires is_container_specification<container_aliases>
  void molecule<container_aliases>::resize(const std::size_t n_atoms, std::size_t n_bonds) {
    mudock::resize(atom_elements, n_atoms);
    mudock::resize(atom_autodock_type, n_atoms);
    mudock::resize(x_coordinates, n_atoms);
    mudock::resize(y_coordinates, n_atoms);
    mudock::resize(z_coordinates, n_atoms);
    mudock::resize(atom_is_aromatic, n_atoms);
    mudock::resize(atom_Rii, n_atoms);
    mudock::resize(atom_vol, n_atoms);
    mudock::resize(atom_solpar, n_atoms);
    mudock::resize(atom_epsii, n_atoms);
    mudock::resize(atom_Rij_hb, n_atoms);
    mudock::resize(atom_epsij_hb, n_atoms);
    mudock::resize(atom_charge, n_atoms);
    mudock::resize(atom_num_hbond, n_atoms);
    mudock::resize(bond_descriptions, n_bonds);
    atoms_size = n_atoms;
    bonds_size = n_bonds;
  }
  template<class container_aliases>
    requires is_container_specification<container_aliases>
  void molecule<container_aliases>::remove_atom(const std::size_t index) {
    // remove the target atom from all the containers
    mudock::remove_atom(atom_elements, index);
    mudock::remove_atom(atom_autodock_type, index);
    mudock::remove_atom(x_coordinates, index);
    mudock::remove_atom(y_coordinates, index);
    mudock::remove_atom(z_coordinates, index);
    mudock::remove_atom(atom_is_aromatic, index);
    mudock::remove_atom(atom_Rii, index);
    mudock::remove_atom(atom_vol, index);
    mudock::remove_atom(atom_solpar, index);
    mudock::remove_atom(atom_epsii, index);
    mudock::remove_atom(atom_Rij_hb, index);
    mudock::remove_atom(atom_epsij_hb, index);
    mudock::remove_atom(atom_charge, index);
    mudock::remove_atom(atom_num_hbond, index);
    atoms_size--;

    // now we need to update the bonds as well
    auto end_loop = std::begin(bond_descriptions) + bonds_size;
    for (auto bond_it{std::begin(bond_descriptions)}; bond_it != end_loop; ++bond_it) {
      auto& source = bond_it->source;
      auto& dest   = bond_it->dest;
      if (source == index || dest == index) {
        std::shift_left(bond_it, end_loop, 1);
        --end_loop;
      } else {
        if (source > index) {
          --source;
        }
        if (dest > index) {
          --dest;
        }
      }
    }
    const auto new_bond_size = std::size_t{end_loop - std::begin(bond_descriptions)};
    mudock::resize(bond_descriptions, new_bond_size);
    bonds_size = new_bond_size;
  }

} // namespace mudock
