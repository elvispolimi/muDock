#include <mudock/chem/bond_types.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<>
  void molecule_type<static_container_type>::resize(const std::size_t n_atoms, std::size_t n_bonds) {
    atoms.resize(n_atoms);
    assert(n_bonds < max_static_bonds());
    num_bonds = static_cast<index_type>(n_bonds);
  }

  template<>
  void molecule_type<dynamic_container_type>::resize(const std::size_t n_atoms, std::size_t n_bonds) {
    atoms.resize(n_atoms);
    bonds.resize(n_bonds, {index_type{0}, index_type{0}, bond_type::SINGLE});
    num_bonds = static_cast<index_type>(n_bonds);
  }

} // namespace mudock
