#pragma once

#include <mudock/batch.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<typename queue_type, typename layout_set = all_scratch_layouts>
  struct docking: public stage<queue_type, layout_set> {
    docking(std::shared_ptr<scratchpad<queue_type, layout_set>> _scratch): stage<queue_type, layout_set>(_scratch) {};
    virtual void prepare(batch<static_molecule>&) = 0;
    virtual void operator()()                     = 0;

    virtual ~docking() = default;
  };
} // namespace mudock
