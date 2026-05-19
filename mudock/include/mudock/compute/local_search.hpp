#pragma once

#include <mudock/batch.hpp>
#include <mudock/molecule.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>

namespace mudock {
  template<typename queue_t, template<typename> typename scoring_t>
    requires std::derived_from<queue_t, queue> && std::derived_from<scoring_t<queue_t>, scoring<queue_t>>
  struct local_search: public stage<queue_t> {
    local_search(std::shared_ptr<scratchpad<queue_t>> _scratch, 
                 std::shared_ptr<scoring_t<queue_t>> _score)
                 : stage<queue_t>(_scratch),
                   score_stage(_score) {};
    virtual void prepare(batch<static_molecule>&) = 0;
    virtual void operator()()                     = 0;

    virtual ~local_search() = default;
    
  protected:
    std::shared_ptr<scoring_t<queue_t>> score_stage;
    size_t iterations;
    size_t convergence_patience;
  };
} // namespace mudock
