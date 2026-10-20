#pragma once

#include <memory>
#include <mudock/batch.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  template<typename queue_type, typename layout_set = all_scratch_layouts>
  struct transform: public stage<queue_type, layout_set> {
    transform(std::shared_ptr<scratchpad<queue_type, layout_set>> _scratch): stage<queue_type, layout_set>(_scratch) {};
    virtual void prepare(batch<static_molecule>&) = 0;
    virtual void operator()()                     = 0;

    virtual ~transform() = default;
  };
} // namespace mudock
