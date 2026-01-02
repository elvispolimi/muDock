#pragma once

#include <functional>
#include <mudock/molecule.hpp>

namespace mudock {

  template<typename molecule_type>
    requires is_molecule<molecule_type>
  void assign_autodock_types(molecule_type&, std::function<void(molecule_type&)> = {});

  template<>
  void assign_autodock_types(static_molecule&, std::function<void(static_molecule&)>);
  template<>
  void assign_autodock_types(dynamic_molecule&, std::function<void(dynamic_molecule&)>);
} // namespace mudock
