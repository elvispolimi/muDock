#include <mudock/chem/bond_types.hpp>
#include <mudock/chem/elements.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<>
  void molecule<static_containers>::resize(const std::size_t n_atoms, std::size_t n_bonds) {
    coordinates.resize(n_atoms);
    assert(n_atoms < max_static_atoms());
    assert(n_bonds < max_static_bonds());
    atom_size  = static_cast<index_type>(n_atoms);
    bonds_size = static_cast<index_type>(n_bonds);
  }

  template<>
  void molecule<dynamic_containers>::resize(const std::size_t n_atoms, std::size_t n_bonds) {
    coordinates.resize(n_atoms);
    atom_elements.resize(n_atoms, element::H);
    bond_descriptions.resize(n_bonds);
    atom_size  = static_cast<index_type>(n_atoms);
    bonds_size = static_cast<index_type>(n_bonds);
  }

} // namespace mudock
