#pragma once

#include <cassert>
#include <memory>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  // define how all the kernel implementations should behave. The idea is that they take the protein as input
  // when the class is constructed, so it can be reused on all the subsequent virtual screening calls. The
  // latter takes one or more ligands in input and provides in output their score
  class kernel_interface {
  protected:
    std::shared_ptr<dynamic_molecule> target_protein;

  public:
    inline kernel_interface(std::shared_ptr<dynamic_molecule>& protein): target_protein(protein) {}
    virtual ~kernel_interface() {}

    virtual void virtual_screen(const std::span<std::unique_ptr<static_molecule>> ligands) = 0;
  };

  // this is a simple wrapper to hide implementation details about the kernel and the remove the burden of
  // dealing with virtual functions
  class kernel {
    std::unique_ptr<kernel_interface> implementation;

  public:
    template<class kernel_type, class... T>
    inline kernel(T&&... args): implementation(std::make_unique<kernel_type>(std::forward(args...))) {}
    inline kernel(): implementation(nullptr) {}
    inline kernel(kernel&&)                 = default;
    inline kernel(const kernel&)            = delete;
    inline ~kernel()                        = default;
    inline kernel& operator=(const kernel&) = delete;
    inline kernel& operator=(kernel&&)      = default;

    inline void virtual_screen(const std::span<std::unique_ptr<static_molecule>> ligands) {
      assert(implementation);
      implementation->virtual_screen(ligands);
    }
  };

} // namespace mudock
