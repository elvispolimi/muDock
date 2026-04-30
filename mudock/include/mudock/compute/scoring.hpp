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
    
    virtual ~scoring() = default;
  
  };

  template<typename queue_type>
  struct differentiable_scoring: public scoring<queue_type> {
    differentiable_scoring(std::shared_ptr<scratchpad<queue_type>> _scratch): scoring<queue_type>(_scratch) {};
    virtual void compute_gradient()               = 0;
    // Using the same structure as operator() for the gradient (using a kernel) instead of having a simple function
    // virtual const gradient& compute_gradient() {
    //   throw std::runtime_error("Gradient not implemented.");
    //   return this->grad;
    // };
  };
} // namespace mudock
