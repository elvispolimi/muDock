#pragma once

#include <memory>
#include <mudock/cuda_implementation/cuda_object.cuh>
#include <mudock/cuda_implementation/cuda_wrapper.cuh>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <random>
#include <span>

namespace mudock {

  class virtual_screen_cuda {
    // this will be the scratchpad memory for the CUDA implementation
    // NOTE: we will implement a random scoring function to test the infrastructure
    std::mt19937 generator;
    std::uniform_real_distribution<fp_type> dist;

  public:
    virtual_screen_cuda(std::shared_ptr<dynamic_molecule>& protein);

    void operator()(std::span<std::unique_ptr<static_molecule>> ligands);
  };
} // namespace mudock
