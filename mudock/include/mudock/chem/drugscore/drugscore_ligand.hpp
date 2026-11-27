#pragma once

#include <functional>
#include <mudock/chem/drugscore/drugscore_layer.hpp>

namespace mudock {
  struct drugscore_ligand: public drugscore_static_layer {
    drugscore_ligand(static_molecule& _molecule, std::function<void(drugscore_static_layer&)> f = {})
        : drugscore_static_layer(_molecule, f) {};
  };
} // namespace mudock
