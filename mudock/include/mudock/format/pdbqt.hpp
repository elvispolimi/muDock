#pragma once

#include <mudock/format/ob_wrapper.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  struct pdbqt {
    static constexpr auto PDBQT_ATOM_TOKEN    = "ATOM";
    static constexpr auto PDBQT_HETATOM_TOKEN = "HETATOM";
  };
} // namespace mudock
