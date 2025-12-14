#pragma once

#include <memory>
#include <mudock/batch.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<typename queue_type>
  struct stage {
    stage(std::shared_ptr<scratchpad<queue_type>> _scratch): scratch(_scratch) {};
    virtual void prepare(batch<static_molecule>&) = 0;
    virtual void operator()()                     = 0;
    void teardown(batch<static_molecule>& b) {
      teardown_impl(b);
      invalid_scratch();
    };

    virtual ~stage() = default;

    void invalid_scratch() { scratch->invalidate(); }

  protected:
    std::shared_ptr<scratchpad<queue_type>> scratch;

    virtual void teardown_impl(batch<static_molecule>&) = 0;
  };
} // namespace mudock
