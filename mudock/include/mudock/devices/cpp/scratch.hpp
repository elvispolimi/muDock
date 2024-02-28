#pragma once

#include <algorithm>
#include <cassert>
#include <mudock/devices/language_types.hpp>
#include <mudock/devices/scratch.hpp>
#include <mudock/molecule.hpp>
#include <span>
#include <vector>

namespace mudock {
  // This is the full specialization of the scratchpad for the CPP language
  template<>
  class language_scratchpad_impl<language_types::CPP>
      : protected language_scratchpad_interface<language_types::CPP> {
  public:
    language_scratchpad_impl(const device_context<language_types::CPP>& d_c)
        : language_scratchpad_interface(d_c){};
    void copy_in(const std::span<std::unique_ptr<static_molecule>> molecules) override {
      this->ligands = molecules;
    };
    void copy_out(std::span<std::unique_ptr<static_molecule>> molecules) override { molecules = ligands; };
    std::span<std::unique_ptr<static_molecule>> get_ligands() { return ligands; };

  private:
    std::span<std::unique_ptr<static_molecule>> ligands;
  };
} // namespace mudock
