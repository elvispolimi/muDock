#pragma once

#include <mudock/chem/autodock_babel_types.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/grid/mdspan.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  template<class container_aliases>
    requires is_container_specification<container_aliases>
  struct autodock_molecule: public molecule<container_aliases> {
    [[nodiscard]] inline auto get_autodock_type() {
      return make_span(atom_autodock_type, molecule<container_aliases>::num_atoms());
    }
    [[nodiscard]] inline auto get_is_aromatic() const {
      return make_span(atom_is_aromatic, molecule<container_aliases>::num_atoms());
    }
    [[nodiscard]] inline auto get_Rii() const {
      return make_span(atom_Rii, molecule<container_aliases>::num_atoms());
    }
    [[nodiscard]] inline auto get_vol() const {
      return make_span(atom_vol, molecule<container_aliases>::num_atoms());
    }
    [[nodiscard]] inline auto get_solpar() const {
      return make_span(atom_solpar, molecule<container_aliases>::num_atoms());
    }
    [[nodiscard]] inline auto get_epsii() const {
      return make_span(atom_epsii, molecule<container_aliases>::num_atoms());
    }
    [[nodiscard]] inline auto get_Rij_hb() const {
      return make_span(atom_Rij_hb, molecule<container_aliases>::num_atoms());
    }
    [[nodiscard]] inline auto get_epsij_hb() const {
      return make_span(atom_epsij_hb, molecule<container_aliases>::num_atoms());
    }
    [[nodiscard]] inline auto get_charge() const {
      return make_span(atom_charge, molecule<container_aliases>::num_atoms());
    }
    [[nodiscard]] inline auto get_num_hbond() const {
      return make_span(atom_num_hbond, molecule<container_aliases>::num_atoms());
    }

    [[nodiscard]] inline auto& autodock_type(const int index) { return atom_autodock_type[index]; }
    [[nodiscard]] inline auto& is_aromatic(const int index) { return atom_is_aromatic[index]; }
    [[nodiscard]] inline auto& Rii(const int index) { return atom_Rii[index]; }
    [[nodiscard]] inline auto& vol(const int index) { return atom_vol[index]; }
    [[nodiscard]] inline auto& solpar(const int index) { return atom_solpar[index]; }
    [[nodiscard]] inline auto& epsii(const int index) { return atom_epsii[index]; }
    [[nodiscard]] inline auto& Rij_hb(const int index) { return atom_Rij_hb[index]; }
    [[nodiscard]] inline auto& epsij_hb(const int index) { return atom_epsij_hb[index]; }
    [[nodiscard]] inline auto& charge(const int index) { return atom_charge[index]; }
    [[nodiscard]] inline auto& num_hbond(const int index) { return atom_num_hbond[index]; }

    [[nodiscard]] inline auto* autodock_type() { return atom_autodock_type.data(); }
    [[nodiscard]] inline auto* is_aromatic() { return atom_is_aromatic.data(); }
    [[nodiscard]] inline auto* Rii() { return atom_Rii.data(); }
    [[nodiscard]] inline auto* vol() { return atom_vol.data(); }
    [[nodiscard]] inline auto* solpar() { return atom_solpar.data(); }
    [[nodiscard]] inline auto* epsii() { return atom_epsii.data(); }
    [[nodiscard]] inline auto* Rij_hb() { return atom_Rij_hb.data(); }
    [[nodiscard]] inline auto* epsij_hb() { return atom_epsij_hb.data(); }
    [[nodiscard]] inline auto* charge() { return atom_charge.data(); }
    [[nodiscard]] inline auto* num_hbond() { return atom_num_hbond.data(); }

    [[nodiscard]] inline const auto& autodock_type(const int index) const {
      return atom_autodock_type[index];
    }
    [[nodiscard]] inline const auto& is_aromatic(const int index) const { return atom_is_aromatic[index]; }
    [[nodiscard]] inline const auto& Rii(const int index) const { return atom_Rii[index]; }
    [[nodiscard]] inline const auto& vol(const int index) const { return atom_vol[index]; }
    [[nodiscard]] inline const auto& solpar(const int index) const { return atom_solpar[index]; }
    [[nodiscard]] inline const auto& epsii(const int index) const { return atom_epsii[index]; }
    [[nodiscard]] inline const auto& Rij_hb(const int index) const { return atom_Rij_hb[index]; }
    [[nodiscard]] inline const auto& epsij_hb(const int index) const { return atom_epsij_hb[index]; }
    [[nodiscard]] inline auto& charge(const int index) const { return atom_charge[index]; }
    [[nodiscard]] inline const auto& num_hbond(const int index) const { return atom_num_hbond[index]; }

    void resize(const int n_atoms, int n_bonds) {
      mudock::resize(atom_autodock_type, n_atoms);
      mudock::resize(atom_is_aromatic, n_atoms);
      mudock::resize(atom_Rii, n_atoms);
      mudock::resize(atom_vol, n_atoms);
      mudock::resize(atom_solpar, n_atoms);
      mudock::resize(atom_epsii, n_atoms);
      mudock::resize(atom_Rij_hb, n_atoms);
      mudock::resize(atom_epsij_hb, n_atoms);
      mudock::resize(atom_charge, n_atoms);
      mudock::resize(atom_num_hbond, n_atoms);
      molecule<container_aliases>::resize(n_atoms, n_bonds);
    }

    void remove_atom(const int index) {
      mudock::remove_atom(atom_autodock_type, index);
      mudock::remove_atom(atom_is_aromatic, index);
      mudock::remove_atom(atom_Rii, index);
      mudock::remove_atom(atom_vol, index);
      mudock::remove_atom(atom_solpar, index);
      mudock::remove_atom(atom_epsii, index);
      mudock::remove_atom(atom_Rij_hb, index);
      mudock::remove_atom(atom_epsij_hb, index);
      mudock::remove_atom(atom_charge, index);
      mudock::remove_atom(atom_num_hbond, index);
      molecule<container_aliases>::remove_atom(index);
    }

    void prepare(std::function<void(autodock_molecule<container_aliases>&)> f = {}) {
      assign_autodock_types(f);
      // fill the atom properties using the autodock force field
      for (int index{0}; index < this->num_atoms(); ++index) {
        const auto& ff_entry = get_description(atom_autodock_type[index]);
        // autodock_type(index) = ff_entry.value;
        Rii(index)       = ff_entry.Rii;
        epsii(index)     = ff_entry.epsii * autodock_parameters::coeff_vdW;
        vol(index)       = ff_entry.vol;
        solpar(index)    = ff_entry.solpar;
        Rij_hb(index)    = ff_entry.Rij_hb;
        epsij_hb(index)  = ff_entry.epsij_hb * autodock_parameters::coeff_hbond;
        num_hbond(index) = ff_entry.hbond;
      }
    }

  private:
    molecule<container_aliases>::template atoms_array_type<autodock_ff> atom_autodock_type;
    molecule<container_aliases>::template atoms_array_type<int>
        atom_is_aromatic; // 1-> aromatic, 0 -> no aromatic
    molecule<container_aliases>::template atoms_array_type<fp_type> atom_Rii;
    molecule<container_aliases>::template atoms_array_type<fp_type> atom_vol;
    molecule<container_aliases>::template atoms_array_type<fp_type> atom_solpar;
    molecule<container_aliases>::template atoms_array_type<fp_type> atom_epsii;
    molecule<container_aliases>::template atoms_array_type<fp_type> atom_Rij_hb;
    molecule<container_aliases>::template atoms_array_type<fp_type> atom_epsij_hb;
    molecule<container_aliases>::template atoms_array_type<fp_type> atom_charge;
    molecule<container_aliases>::template atoms_array_type<int> atom_num_hbond;
    void assign_autodock_types(std::function<void(autodock_molecule<container_aliases>&)> f = {});
  };

  using autodock_dynamic_molecule = autodock_molecule<dynamic_containers>;
  using autodock_static_molecule  = autodock_molecule<static_containers>;

  template<class T>
  concept is_autodock_molecule = (std::same_as<std::remove_cvref_t<T>, autodock_static_molecule> ||
                                  std::same_as<std::remove_cvref_t<T>, autodock_dynamic_molecule>);
  template<class T>
  concept derived_from_autodock_molecule =
      (std::derived_from<T, autodock_static_molecule> || std::derived_from<T, autodock_dynamic_molecule>);
} // namespace mudock
