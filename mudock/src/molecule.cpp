#include <mudock/chem/bond_types.hpp>
#include <mudock/chem/elements.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<>
  void molecule<static_container_type>::resize(const std::size_t n_atoms, std::size_t n_bonds) {
    coordinates.resize(n_atoms);
    assert(n_atoms < max_static_atoms());
    assert(n_bonds < max_static_bonds());
    num_atoms = static_cast<index_type>(n_atoms);
    num_bonds = static_cast<index_type>(n_bonds);
  }

  template<>
  void molecule<dynamic_container_type>::resize(const std::size_t n_atoms, std::size_t n_bonds) {
    coordinates.resize(n_atoms);
    elements.resize(n_atoms, element::H);
    bonds.resize(n_bonds);
    num_atoms = static_cast<index_type>(n_atoms);
    num_bonds = static_cast<index_type>(n_bonds);
  }

} // namespace mudock
