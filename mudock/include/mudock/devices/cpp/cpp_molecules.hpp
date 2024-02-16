#pragma once

#include <algorithm>
#include <cassert>
#include <mudock/devices/language_types.hpp>
#include <mudock/devices/molecules.hpp>
#include <mudock/molecule.hpp>
#include <span>
#include <vector>

namespace mudock {
  template<>
  class molecules_scratchpad<language_types::CPP> {
  public:
    molecules_scratchpad(){};
    inline void copy_in(const std::span<std::unique_ptr<static_molecule>>& molecules,
                        const device_context<language_types::CPP>&) {
      this->ligands = molecules;
    };
    inline void copy_out(std::span<std::unique_ptr<static_molecule>> molecules,
                         const device_context<language_types::CPP>&) {
      molecules = ligands;
    };
    std::span<std::unique_ptr<static_molecule>> get_ligands() { return ligands; };

  private:
    std::span<std::unique_ptr<static_molecule>> ligands;
  };
} // namespace mudock
