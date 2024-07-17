#pragma once

#include <memory>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <random>
#include <span>
#include <string>

namespace mudock {

  class virtual_screen_cpp {
    // this will be the scratchpad memory for the cpp implementation
    // NOTE: we will implement a random scoring function to test the infrastructure
    std::mt19937 generator;
    std::uniform_real_distribution<fp_type> dist;

  public:
    virtual_screen_cpp(std::shared_ptr<dynamic_molecule>& protein);

    void operator()(static_molecule& ligand);
  };

} // namespace mudock
