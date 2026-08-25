#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/chem/x_score_xtool_types.hpp>
#include <mudock/chem/x_score_xlogp_types.hpp>
#include <mudock/chem/residue.hpp>
#include <mudock/chem/elements.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/residue_types.hpp>
#include <mudock/chem/sybyl_atom_types.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/molecule/bond.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/molecule/graph.hpp>
#include <mudock/molecule/properties.hpp>
#include <mudock/molecule/property_table.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {
  static constexpr std::string_view DEFAULT_RESIDUE_NAME = "LIG";
  static constexpr int DEFAULT_RESIDUE_ID                = 1;

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
    atoms_array_type<fp_type> x_coordinates;
    atoms_array_type<fp_type> y_coordinates;
    atoms_array_type<fp_type> z_coordinates;
    int atoms_size = int{0};

    // the intra-molecular connections
    bonds_array_type<bond> bond_descriptions;
    int bonds_size = int{0};

    atoms_array_type<std::string> atom_names;
    atoms_array_type<sybyl_atom_type> atom_sybyl_types;
    atoms_array_type<int> residue_ids;
    atoms_array_type<std::string> residue_names;
    atoms_array_type<residue_type> atom_residue_types;

    atoms_array_type<autodock_ff> atom_autodock_type;

    // protein-specific fields (x_score layer)
    // renamed from atom_residue_types to avoid collision with the residue_type enum field above
    atoms_array_type<residue> atom_residue;

    // ligand-specific fields
    atoms_array_type<xtool_ff> atom_types;

  

    atoms_array_type<int> atom_is_aromatic;
    atoms_array_type<fp_type> atom_charge;
    atoms_array_type<int> atom_num_hbond;

  public:
    // functions to manage the geometry of a molecule
    void resize(const int n_atoms, int n_bonds);
    void remove_atom(const int index);
    [[nodiscard]] constexpr auto num_atoms() const { return atoms_size; }
    [[nodiscard]] constexpr auto num_bonds() const { return bonds_size; }
    constexpr auto num_rotamers() const {
      return std::count_if(std::begin(bond_descriptions), std::end(bond_descriptions), [](const bond& b) {
        return b.can_rotate;
      });
    }

    point<fp_type, 3> get_center() const {
      // FIX ME I don't like it here, I do not undestand why it is needed by ADT
      // find the protein boundaries and its center
      auto min = point3D{std::ranges::min(x_coordinates),
                         std::ranges::min(y_coordinates),
                         std::ranges::min(z_coordinates)};
      std::transform(min.begin(), min.end(), min.begin(), [](const auto p) {
        return std::floor((p - (grid_spacing + cutoff_distance)) * inv_spacing) / inv_spacing;
      });
      auto max = point3D{std::ranges::max(x_coordinates),
                         std::ranges::max(y_coordinates),
                         std::ranges::max(z_coordinates)};
      std::transform(max.begin(), max.end(), max.begin(), [](const auto p) {
        return std::ceil((p + (grid_spacing + cutoff_distance)) * inv_spacing) / inv_spacing;
      });
      return {max.difference(min).divide({fp_type{2}}).add(min)};
    }

    // void prepare() {};

    // a container that we can use to store key-value properties, e.g. its name
    property_map properties;

    // utility functions to get information about the bonds
    [[nodiscard]] inline auto get_bonds() { return std::span(std::begin(bond_descriptions), bonds_size); }
    [[nodiscard]] inline auto get_bonds() const {
      return std::span(std::cbegin(bond_descriptions), bonds_size);
    }
    [[nodiscard]] inline auto& bonds(int index) { return bond_descriptions[index]; }
    [[nodiscard]] inline const auto& bonds(int index) const { return bond_descriptions[index]; }

    // utility functions to get the span of the whole molecule (read + write)
    [[nodiscard]] inline auto get_autodock_type() { return make_span(atom_autodock_type, atoms_size); }
    [[nodiscard]] inline auto get_residue_types() { return make_span(atom_residue, atoms_size); }
    [[nodiscard]] inline auto get_atom_type() { return make_span(atom_types, atoms_size); }
    [[nodiscard]] inline auto get_elements() { return make_span(atom_elements, atoms_size); }
    [[nodiscard]] inline auto get_x() { return make_span(x_coordinates, atoms_size); }
    [[nodiscard]] inline auto get_y() { return make_span(y_coordinates, atoms_size); }
    [[nodiscard]] inline auto get_z() { return make_span(z_coordinates, atoms_size); }
    [[nodiscard]] inline auto get_autodock_type() const { return make_span(atom_autodock_type, atoms_size); }
    [[nodiscard]] inline auto get_is_aromatic() const { return make_span(atom_is_aromatic, atoms_size); }
    [[nodiscard]] inline auto get_charge() const { return make_span(atom_charge, atoms_size); }
    [[nodiscard]] inline auto get_num_hbond() const { return make_span(atom_num_hbond, atoms_size); }
    [[nodiscard]] inline auto get_atom_residue_types() const { return make_span(atom_residue_types, atoms_size); }
    [[nodiscard]] inline auto get_atom_name() const { return make_span(atom_names, atoms_size); }
    [[nodiscard]] inline auto get_atom_type() const { return make_span(atom_types, atoms_size); }

    // utility functions to get the span of the whole molecule (read only)
    [[nodiscard]] inline auto get_elements() const { return make_span(atom_elements, atoms_size); }
    // [[nodiscard]] inline auto get_autodock_type() const { return make_span(atom_autodock_type, atoms_size); }
    [[nodiscard]] inline auto get_x() const { return make_span(x_coordinates, atoms_size); }
    [[nodiscard]] inline auto get_y() const { return make_span(y_coordinates, atoms_size); }
    [[nodiscard]] inline auto get_z() const { return make_span(z_coordinates, atoms_size); }

    // utility functions to get the ref to an atom element (read + write)
    [[nodiscard]] inline auto& autodock_type(const int index) { return atom_autodock_type[index]; }
    [[nodiscard]] inline auto& is_aromatic(const int index) { return atom_is_aromatic[index]; }
    [[nodiscard]] inline auto& elements(const int index) { return atom_elements[index]; }
    [[nodiscard]] inline auto& x(const int index) { return x_coordinates[index]; }
    [[nodiscard]] inline auto& y(const int index) { return y_coordinates[index]; }
    [[nodiscard]] inline auto& z(const int index) { return z_coordinates[index]; }
    [[nodiscard]] inline auto& charge(const int index) { return atom_charge[index]; }
    [[nodiscard]] inline auto& num_hbond(const int index) { return atom_num_hbond[index]; }
    [[nodiscard]] inline auto& residue_types(const int index) { return atom_residue[index]; }
    [[nodiscard]] inline auto& atom_name(const int index) { return atom_names[index]; }
    [[nodiscard]] inline auto& atom_type(const int index) { return atom_types[index]; }

    // utility functions to get the ref to an atom element (read + write)
    [[nodiscard]] inline auto* autodock_type() { return atom_autodock_type.data(); }
    [[nodiscard]] inline auto* is_aromatic() { return atom_is_aromatic.data(); }
    [[nodiscard]] inline auto* elements() { return atom_elements.data(); }
    [[nodiscard]] inline auto* x() { return x_coordinates.data(); }
    [[nodiscard]] inline auto* y() { return y_coordinates.data(); }
    [[nodiscard]] inline auto* z() { return z_coordinates.data(); }
    [[nodiscard]] inline auto* charge() { return atom_charge.data(); }
    [[nodiscard]] inline auto* num_hbond() { return atom_num_hbond.data(); }
    [[nodiscard]] inline auto* residue_types() { return atom_residue.data(); }
    [[nodiscard]] inline auto* atom_name() { return atom_names.data(); }
    [[nodiscard]] inline auto* atom_type() { return atom_types.data(); }

    // utility functions to get the ref to an atom element (read only)
    [[nodiscard]] inline const auto& autodock_type(const int index) const {
      return atom_autodock_type[index];
    }
    [[nodiscard]] inline const auto& is_aromatic(const int index) const { return atom_is_aromatic[index]; }
    [[nodiscard]] inline const auto& elements(const int index) const { return atom_elements[index]; }
    [[nodiscard]] inline const auto& x(const int index) const { return x_coordinates[index]; }
    [[nodiscard]] inline const auto& y(const int index) const { return y_coordinates[index]; }
    [[nodiscard]] inline const auto& z(const int index) const { return z_coordinates[index]; }
    [[nodiscard]] inline const auto& charge(const int index) const { return atom_charge[index]; }
    [[nodiscard]] inline const auto& num_hbond(const int index) const { return atom_num_hbond[index]; }
    [[nodiscard]] inline const auto& residue_types(const int index) const { return atom_residue[index]; }
    [[nodiscard]] inline const auto& atom_name(const int index) const { return atom_names[index]; }
    [[nodiscard]] inline const auto& atom_type(const int index) const { return atom_types[index]; }

    [[nodiscard]] inline auto& sybyl_type(const int index) { return atom_sybyl_types[index]; }
    [[nodiscard]] inline const auto& sybyl_type(const int index) const { return atom_sybyl_types[index]; }

    [[nodiscard]] inline auto& residue_id(const int index) { return residue_ids[index]; }
    [[nodiscard]] inline const auto& residue_id(const int index) const { return residue_ids[index]; }

    [[nodiscard]] inline auto& residue_name(const int index) { return residue_names[index]; }
    [[nodiscard]] inline const auto& residue_name(const int index) const { return residue_names[index]; }

    [[nodiscard]] inline auto& atom_residue_type(const int index) { return atom_residue_types[index]; }
    [[nodiscard]] inline const auto& atom_residue_type(const int index) const {
      return atom_residue_types[index];
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
  template<class T>
  concept derived_from_molecule =
      (std::derived_from<T, static_molecule> || std::derived_from<T, dynamic_molecule>);

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<class container_aliases>
    requires is_container_specification<container_aliases>
  void molecule<container_aliases>::resize(const int n_atoms, int n_bonds) {
    mudock::resize(atom_elements, n_atoms);
    mudock::resize(x_coordinates, n_atoms);
    mudock::resize(y_coordinates, n_atoms);
    mudock::resize(z_coordinates, n_atoms);
    mudock::resize(bond_descriptions, n_bonds);
    mudock::resize(atom_autodock_type, n_atoms);
    mudock::resize(atom_residue, n_atoms);
    mudock::resize(atom_types, n_atoms);
    mudock::resize(atom_is_aromatic, n_atoms);
    mudock::resize(atom_charge, n_atoms);
    mudock::resize(atom_num_hbond, n_atoms);
    mudock::resize(atom_names, n_atoms);
    mudock::resize(atom_sybyl_types, n_atoms);

    mudock::resize(residue_ids, n_atoms);
    mudock::resize(residue_names, n_atoms);
    mudock::resize(atom_residue_types, n_atoms);

    std::fill(std::begin(atom_sybyl_types), std::begin(atom_sybyl_types) + n_atoms, sybyl_atom_type::UNKNOWN);

    std::fill(std::begin(residue_ids), std::begin(residue_ids) + n_atoms, DEFAULT_RESIDUE_ID);
    std::fill(std::begin(residue_names), std::begin(residue_names) + n_atoms, DEFAULT_RESIDUE_NAME);

    std::fill(std::begin(atom_residue_types), std::begin(atom_residue_types) + n_atoms, residue_type::LIG);
    atoms_size = n_atoms;
    bonds_size = n_bonds;
  }
  template<class container_aliases>
    requires is_container_specification<container_aliases>
  void molecule<container_aliases>::remove_atom(const int index) {
    // remove the target atom from all the containers
    mudock::remove_atom(atom_elements, index);
    mudock::remove_atom(x_coordinates, index);
    mudock::remove_atom(y_coordinates, index);
    mudock::remove_atom(z_coordinates, index);
    mudock::remove_atom(atom_residue, index);
    mudock::remove_atom(atom_types, index);
    mudock::resize(bond_descriptions, index);
    mudock::resize(atom_autodock_type, index);
    mudock::resize(atom_is_aromatic, index);
    mudock::resize(atom_charge, index);
    mudock::resize(atom_num_hbond, index);
    mudock::remove_atom(atom_names, index);
    mudock::remove_atom(atom_sybyl_types, index);

    mudock::remove_atom(residue_ids, index);
    mudock::remove_atom(residue_names, index);
    mudock::remove_atom(atom_residue_types, index);
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
    const auto new_bond_size = int{end_loop - std::begin(bond_descriptions)};
    mudock::resize(bond_descriptions, new_bond_size);
    bonds_size = new_bond_size;
  }

} // namespace mudock
