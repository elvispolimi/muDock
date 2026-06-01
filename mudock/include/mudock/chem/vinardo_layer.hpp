#pragma once

#include <cstdint>

#include <mudock/chem/assign_vinardo_type.hpp>
#include <mudock/chem/molecule_layer.hpp>
#include <mudock/chem/vinardo_type.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>

/*
This layer stores the atom-level static properties needed to build a Vinardo
view of a molecule and to support later preprocessing and scoring.
*/

namespace mudock {
  template<class container_aliases>
    requires is_container_specification<container_aliases>
  struct vinardo_layer: public molecule_layer<container_aliases> {

    template<typename T>
    using atoms_array_type = container_aliases::template atoms_size<T>;
    using molecule_type    = molecule<container_aliases>;

    vinardo_layer(molecule_type& molecule): molecule_layer<container_aliases>(molecule) {
      const auto num_atoms = molecule.num_atoms();
      mudock::resize(atom_vinardo_type, num_atoms);
      mudock::resize(atom_radius, num_atoms);
      mudock::resize(atom_is_hydrophobic, num_atoms);
      mudock::resize(atom_is_hbond_donor, num_atoms);
      mudock::resize(atom_is_hbond_acceptor, num_atoms);
      prepare(molecule);
    }

    [[nodiscard]] inline auto get_vinardo_type() const {
      return make_span(atom_vinardo_type, this->get_base_molecule().num_atoms());
    }

    [[nodiscard]] inline auto get_is_hydrophobic() const {
      return make_span(atom_is_hydrophobic, this->get_base_molecule().num_atoms());
    }

    [[nodiscard]] inline auto get_is_hbond_donor() const {
      return make_span(atom_is_hbond_donor, this->get_base_molecule().num_atoms());
    }

    [[nodiscard]] inline auto get_is_hbond_acceptor() const {
      return make_span(atom_is_hbond_acceptor, this->get_base_molecule().num_atoms());
    }

    [[nodiscard]] inline auto get_radius() const {
      return make_span(atom_radius, this->get_base_molecule().num_atoms());
    }

    [[nodiscard]] inline auto num_atoms() const { return this->get_base_molecule().num_atoms(); }

  private:
    void prepare(molecule_type& molecule) {
      assign_vinardo_type(molecule,
                          this->atom_vinardo_type,
                          this->atom_radius,
                          this->atom_is_hydrophobic,
                          this->atom_is_hbond_donor,
                          this->atom_is_hbond_acceptor);
    }

    atoms_array_type<vinardo_atom_type> atom_vinardo_type;
    atoms_array_type<fp_type> atom_radius;

    atoms_array_type<std::uint8_t> atom_is_hydrophobic;
    atoms_array_type<std::uint8_t> atom_is_hbond_donor;
    atoms_array_type<std::uint8_t> atom_is_hbond_acceptor;
  };

} // namespace mudock
