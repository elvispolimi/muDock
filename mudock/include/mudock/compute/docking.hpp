#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<typename queue_type>
  struct docking: public stage<queue_type> {
    docking(std::shared_ptr<scratchpad<queue_type>> _scratch): stage<queue_type>(_scratch) {};
    virtual void prepare(batch<static_molecule>&) = 0;
    virtual void operator()()                     = 0;

    virtual ~docking() = default;
  };
} // namespace mudock
