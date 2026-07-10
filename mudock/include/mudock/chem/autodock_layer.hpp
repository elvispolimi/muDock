#pragma once

#include <mudock/chem/assign_autodock_types.hpp>
#include <mudock/chem/autodock_babel_types.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/molecule_layer.hpp>
#include <mudock/grid/mdspan.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  template<class container_aliases>
    requires is_container_specification<container_aliases>
  struct autodock_layer: public molecule_layer<container_aliases> {
    template<typename T>
    using atoms_array_type = container_aliases::template atoms_size<T>;
    template<typename T>
    using bonds_array_type = container_aliases::template bonds_size<T>;

    autodock_layer(molecule<container_aliases>& _molecule,
                   std::function<void(molecule<container_aliases>&)> f = {})
        : molecule_layer<container_aliases>(_molecule) {
      const auto num_atoms = _molecule.num_atoms();
      mudock::resize(atom_Rii, num_atoms);
      mudock::resize(atom_vol, num_atoms);
      mudock::resize(atom_solpar, num_atoms);
      mudock::resize(atom_epsii, num_atoms);
      mudock::resize(atom_Rij_hb, num_atoms);
      mudock::resize(atom_epsij_hb, num_atoms);
      prepare(f);
    };
    autodock_layer(molecule<container_aliases>& _molecule,
                   std::function<void(autodock_layer<container_aliases>&)> f)
        : molecule_layer<container_aliases>(_molecule) {
      const auto num_atoms = _molecule.num_atoms();
      mudock::resize(atom_Rii, num_atoms);
      mudock::resize(atom_vol, num_atoms);
      mudock::resize(atom_solpar, num_atoms);
      mudock::resize(atom_epsii, num_atoms);
      mudock::resize(atom_Rij_hb, num_atoms);
      mudock::resize(atom_epsij_hb, num_atoms);
      f(*this);
    };

    [[nodiscard]] inline auto get_Rii() const {
      return make_span(atom_Rii, *this->get_base_molecule().num_atoms());
    }
    [[nodiscard]] inline auto get_vol() const {
      return make_span(atom_vol, this->get_base_molecule().num_atoms());
    }
    [[nodiscard]] inline auto get_solpar() const {
      return make_span(atom_solpar, this->get_base_molecule().num_atoms());
    }
    [[nodiscard]] inline auto get_epsii() const {
      return make_span(atom_epsii, this->get_base_molecule().num_atoms());
    }
    [[nodiscard]] inline auto get_Rij_hb() const {
      return make_span(atom_Rij_hb, this->get_base_molecule().num_atoms());
    }
    [[nodiscard]] inline auto get_epsij_hb() const {
      return make_span(atom_epsij_hb, this->get_base_molecule().num_atoms());
    }

    [[nodiscard]] inline auto& Rii(const int index) { return atom_Rii[index]; }
    [[nodiscard]] inline auto& vol(const int index) { return atom_vol[index]; }
    [[nodiscard]] inline auto& solpar(const int index) { return atom_solpar[index]; }
    [[nodiscard]] inline auto& epsii(const int index) { return atom_epsii[index]; }
    [[nodiscard]] inline auto& Rij_hb(const int index) { return atom_Rij_hb[index]; }
    [[nodiscard]] inline auto& epsij_hb(const int index) { return atom_epsij_hb[index]; }

    [[nodiscard]] inline auto* Rii() { return atom_Rii.data(); }
    [[nodiscard]] inline auto* vol() { return atom_vol.data(); }
    [[nodiscard]] inline auto* solpar() { return atom_solpar.data(); }
    [[nodiscard]] inline auto* epsii() { return atom_epsii.data(); }
    [[nodiscard]] inline auto* Rij_hb() { return atom_Rij_hb.data(); }
    [[nodiscard]] inline auto* epsij_hb() { return atom_epsij_hb.data(); }

    [[nodiscard]] inline const auto& Rii(const int index) const { return atom_Rii[index]; }
    [[nodiscard]] inline const auto& vol(const int index) const { return atom_vol[index]; }
    [[nodiscard]] inline const auto& solpar(const int index) const { return atom_solpar[index]; }
    [[nodiscard]] inline const auto& epsii(const int index) const { return atom_epsii[index]; }
    [[nodiscard]] inline const auto& Rij_hb(const int index) const { return atom_Rij_hb[index]; }
    [[nodiscard]] inline const auto& epsij_hb(const int index) const { return atom_epsij_hb[index]; }

    void resize(const int n_atoms, int n_bonds) {
      mudock::resize(atom_Rii, n_atoms);
      mudock::resize(atom_vol, n_atoms);
      mudock::resize(atom_solpar, n_atoms);
      mudock::resize(atom_epsii, n_atoms);
      mudock::resize(atom_Rij_hb, n_atoms);
      mudock::resize(atom_epsij_hb, n_atoms);
      molecule<container_aliases>::resize(n_atoms, n_bonds);
    }

    void remove_atom(const int index) {
      mudock::remove_atom(atom_Rii, index);
      mudock::remove_atom(atom_vol, index);
      mudock::remove_atom(atom_solpar, index);
      mudock::remove_atom(atom_epsii, index);
      mudock::remove_atom(atom_Rij_hb, index);
      mudock::remove_atom(atom_epsij_hb, index);
      molecule<container_aliases>::remove_atom(index);
    }
    [[nodiscard]] inline auto num_atoms() const { return this->get_base_molecule().num_atoms(); }
    [[nodiscard]] inline auto num_bonds() const { return this->get_base_molecule().num_bonds(); }
    [[nodiscard]] inline auto num_rotamers() const { return this->get_base_molecule().num_rotamers(); }
    [[nodiscard]] inline auto get_bonds() const { return this->get_base_molecule().get_bonds(); }
    [[nodiscard]] inline auto get_elements() const { return this->get_base_molecule().get_elements(); }
    [[nodiscard]] inline auto get_x() const { return this->get_base_molecule().get_x(); }
    [[nodiscard]] inline auto get_y() const { return this->get_base_molecule().get_y(); }
    [[nodiscard]] inline auto get_z() const { return this->get_base_molecule().get_z(); }

    [[nodiscard]] inline auto num_atoms() { return this->get_base_molecule().num_atoms(); }
    [[nodiscard]] inline auto num_bonds() { return this->get_base_molecule().num_bonds(); }
    [[nodiscard]] inline auto num_rotamers() { return this->get_base_molecule().num_rotamers(); }
    [[nodiscard]] inline auto get_bonds() { return this->get_base_molecule().get_bonds(); }
    [[nodiscard]] inline auto get_elements() { return this->get_base_molecule().get_elements(); }
    [[nodiscard]] inline auto get_x() { return this->get_base_molecule().get_x(); }
    [[nodiscard]] inline auto get_y() { return this->get_base_molecule().get_y(); }
    [[nodiscard]] inline auto get_z() { return this->get_base_molecule().get_z(); }

  private:
    atoms_array_type<fp_type> atom_Rii;
    atoms_array_type<fp_type> atom_vol;
    atoms_array_type<fp_type> atom_solpar;
    atoms_array_type<fp_type> atom_epsii;
    atoms_array_type<fp_type> atom_Rij_hb;
    atoms_array_type<fp_type> atom_epsij_hb;

    void prepare(std::function<void(molecule<container_aliases>&)> f = {}) {
      assign_autodock_types((*this)(), f);
      // fill the atom properties using the autodock force field
      for (int index{0}; index < this->get_base_molecule().num_atoms(); ++index) {
        const auto& ff_entry = get_description(this->get_base_molecule().autodock_type(index));
        Rii(index)           = ff_entry.Rii;
        epsii(index)         = ff_entry.epsii * autodock_parameters::coeff_vdW;
        vol(index)           = ff_entry.vol;
        solpar(index)        = ff_entry.solpar;
        Rij_hb(index)        = ff_entry.Rij_hb;
        epsij_hb(index)      = ff_entry.epsij_hb * autodock_parameters::coeff_hbond;
        this->get_base_molecule().num_hbond(index) = ff_entry.hbond;
      }
    }
  };

  using autodock_dynamic_layer = autodock_layer<dynamic_containers>;
  using autodock_static_layer  = autodock_layer<static_containers>;

  template<class T>
  concept is_autodock_layer = (std::same_as<std::remove_cvref_t<T>, autodock_static_layer> ||
                               std::same_as<std::remove_cvref_t<T>, autodock_dynamic_layer>);
  template<class T>
  concept derived_from_autodock_layer =
      (std::derived_from<T, autodock_static_layer> || std::derived_from<T, autodock_dynamic_layer>);
} // namespace mudock
