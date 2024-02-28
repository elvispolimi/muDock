#pragma once

#include <algorithm>
#include <cassert>
#include <mudock/devices/kernel_types.hpp>
#include <mudock/devices/scratch.hpp>
#include <mudock/molecule.hpp>
#include <span>
#include <vector>

namespace mudock {
  // This is the full specialization of the scratchpad for the CPP kernel
  template<>
  class kernel_scratchpad_impl<kernel_type::CPP>: protected kernel_scratchpad_interface<kernel_type::CPP> {
  public:
    kernel_scratchpad_impl(const device_context<kernel_type::CPP>& d_c): kernel_scratchpad_interface(d_c){};
    void copy_in(const std::span<std::unique_ptr<static_molecule>> i_molecules) override {
      this->ligands = i_molecules;
    };
    void copy_out(std::span<std::unique_ptr<static_molecule>> o_molecules) override {
      o_molecules = ligands;
    };
    [[nodiscard]] std::span<std::unique_ptr<static_molecule>> get_ligands() { return ligands; };

  private:
    std::span<std::unique_ptr<static_molecule>> ligands;
  };
} // namespace mudock
