#pragma once

#include <array>
#include <concepts>
#include <mudock/molecule/constraints.hpp>
#include <vector>

namespace mudock {

  // this set of type alias describe a molecule that allocates a fixed size of memory to hold atoms and bonds
  // properties. On the bright side, we can reduce the number of allocs, but we need to introduce constraints
  // on the maximum number of atoms and bonds that we can process. We typically use this setup for ligands.
  struct static_containers {
    template<typename T>
    using atoms_size = std::array<T, max_static_atoms()>;

    template<typename T>
    using bonds_size = std::array<T, max_static_bonds()>;

    template<typename T>
    using fragments_size = std::array<T, max_static_atoms() * max_static_bonds()>;
  };

  // this set of type alas describe a molecule that allocates the required memory at runtime.It is the most
  // flexible approch, but it introduces runtime overheads. We typically use this setup for the proteins.
  struct dynamic_containers {
    template<typename T>
    using atoms_size = std::vector<T>;

    template<typename T>
    using bonds_size = std::vector<T>;

    template<typename T>
    using fragments_size = std::vector<T>;
  };

  // this concepts tell if a type is a specification of the container abstraction used in a class.
  template<class T>
  concept is_container_specification = (std::same_as<std::remove_cvref_t<T>, static_containers> ||
                                        std::same_as<std::remove_cvref_t<T>, dynamic_containers>);

} // namespace mudock
