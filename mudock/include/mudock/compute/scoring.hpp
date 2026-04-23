#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/molecule.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>

namespace mudock {
  template<typename queue_type>
  struct scoring: public stage<queue_type> {
    scoring(std::shared_ptr<scratchpad<queue_type>> _scratch): stage<queue_type>(_scratch) {};
    virtual void prepare(batch<static_molecule>&) = 0;
    virtual void operator()()                     = 0;
    // Like the operator computes the score, this should compute the gradient
    chromosome gradient;
    virtual const chromosome& compute_gradient() {
      throw std::runtime_error("Gradient not implemented.");
      return this->gradient;
    };
    //virtual chromosome compute_gradient(const chromosome& pose, const static_molecule& ligand, const dynamic_molecule& protein) = 0;

    virtual ~scoring() = default;
  };
} // namespace mudock
