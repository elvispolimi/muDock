#pragma once

#include <mudock/batch.hpp>
#include <mudock/molecule.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>

namespace mudock {
    // TODO should it extend stage?
  struct local_search: public stage<queue_type> {
    local_search(std::shared_ptr<scratchpad<queue_type>> _scratch): stage<queue_type>(_scratch) {};
    virtual void prepare(batch<static_molecule>&) = 0;
    virtual void operator()()                     = 0;

    virtual ~local_search() = default;
    //private:
    // scoring_t<queue_t>& score_stage;???
  };
} // namespace mudock
