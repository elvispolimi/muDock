#pragma once

#include <memory>
#include <mudock/batch.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  template<typename queue_type>
  struct convergence: public stage<queue_type> {
    convergence(std::shared_ptr<scratchpad<queue_type>> _scratch): stage<queue_type>(_scratch) {};
    virtual void prepare(batch<static_molecule>&) = 0;
    virtual void operator()()                     = 0;

    virtual ~convergence() = default;
  };
} // namespace mudock
