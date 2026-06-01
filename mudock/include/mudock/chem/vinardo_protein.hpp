#pragma once

#include <mudock/chem/vinardo_layer.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  struct vinardo_protein: public vinardo_layer<dynamic_containers> {
    vinardo_protein(dynamic_molecule& protein): vinardo_layer<dynamic_containers>(protein) {}
  };
} // namespace mudock
