#pragma once

#include <mudock/molecule.hpp>
#include <mudock/molecule/containers.hpp>

namespace mudock {
  template<class container_aliases>
    requires is_container_specification<container_aliases>
  struct molecule_layer {
    template<typename T>
    using atoms_array_type = container_aliases::template atoms_size<T>;
    template<typename T>
    using bonds_array_type = container_aliases::template bonds_size<T>;

    molecule_layer(molecule<container_aliases>& mol): base_molecule(mol) {};

    auto& get_base_molecule() const { return base_molecule; };
    auto& operator()() const { return get_base_molecule(); };

    auto& get_base_molecule() { return base_molecule; };
    auto& operator()() { return get_base_molecule(); };

  private:
    molecule<container_aliases>& base_molecule;
  };
} // namespace mudock
