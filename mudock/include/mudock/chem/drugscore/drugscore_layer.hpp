#pragma once

#include <mudock/chem/drugscore/drugscore_types.hpp>
#include <mudock/chem/molecule_layer.hpp>
#include <mudock/type_alias.hpp>
#include <numeric>

namespace mudock {
  template<class container_aliases>
    requires is_container_specification<container_aliases>
  struct drugscore_layer: public molecule_layer<container_aliases> {
    template<typename T>
    using atoms_array_type = container_aliases::template atoms_size<T>;

    drugscore_layer(molecule<container_aliases>& _molecule,
                    std::function<void(drugscore_layer<container_aliases>&)> f = {})
        : molecule_layer<container_aliases>(_molecule) {
      const auto num_atoms = _molecule.num_atoms();
      mudock::resize(atom_drugscore_type, num_atoms);
      prepare(f);
    };

    [[nodiscard]] inline auto get_drugscore_type() const {
      return make_span(atom_drugscore_type, this->get_base_molecule().num_atoms());
    }

    [[nodiscard]] inline auto& drugscore_type(const int index) { return atom_drugscore_type[index]; }

    [[nodiscard]] inline auto* drugscore_type() { return atom_drugscore_type.data(); }
    [[nodiscard]] inline const auto& drugscore_type(const int index) const {
      return atom_drugscore_type[index];
    }

    void resize(const int n_atoms, int n_bonds) {
      this->get_base_molecule().resize(n_atoms, n_bonds);
      mudock::resize(atom_drugscore_type, n_atoms);
    }

    void remove_atom(const int index) {
      this->get_base_molecule().remove_atom(index);
      mudock::remove_atom(atom_drugscore_type, index);
    }

    void permute(const unsigned int* permutation, const unsigned int n_atoms) {
      this->get_base_molecule().permute(permutation, n_atoms);
      mudock::permute(atom_drugscore_type, permutation, n_atoms);
    }

  private:
    atoms_array_type<ds_type> atom_drugscore_type;

    void assign_drugscore_types(std::function<void(drugscore_layer<container_aliases>&)> f = {});
    void prepare(std::function<void(drugscore_layer<container_aliases>&)> f = {}) {
      assign_drugscore_types(f);
    }
  };

  using drugscore_dynamic_layer = drugscore_layer<dynamic_containers>;
  using drugscore_static_layer  = drugscore_layer<static_containers>;

  template<class T>
  concept is_drugscore_layer = (std::same_as<std::remove_cvref_t<T>, drugscore_dynamic_layer> ||
                                std::same_as<std::remove_cvref_t<T>, drugscore_static_layer>);
  template<class T>
  concept derived_from_drugscore_layer =
      (std::derived_from<T, drugscore_dynamic_layer> || std::derived_from<T, drugscore_static_layer>);

} // namespace mudock
