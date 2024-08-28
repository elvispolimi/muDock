#pragma once

#include <memory>
#include <mudock/batch.hpp>
#include <mudock/cuda_implementation/cuda_object.cuh>
#include <mudock/cuda_implementation/cuda_wrapper.cuh>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <random>

namespace mudock {

  class virtual_screen_cuda {
    // this will be the scratchpad memory for the CUDA implementation
    // NOTE: we will implement a random scoring function to test the infrastructure
    std::mt19937 generator;
    std::uniform_real_distribution<fp_type> dist;

    // the configuration of the GA algorithm
    knobs configuration;

  public:
    virtual_screen_cuda(std::shared_ptr<dynamic_molecule>& protein, const knobs k);

    void operator()(batch& incoming_batch);
  };
} // namespace mudock
