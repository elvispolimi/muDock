#pragma once

#include <mudock/batch.hpp>
#include <mudock/data_formats/data_store.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<typename queue_type, typename layout_set = scoring_data_layouts>
  struct scoring: public stage<queue_type, layout_set> {
    scoring(std::shared_ptr<data_store<queue_type, layout_set>> _data): stage<queue_type, layout_set>(_data) {};
    virtual void prepare(batch<static_molecule>&) = 0;
    virtual void operator()()                     = 0;

    virtual ~scoring() = default;
  };
} // namespace mudock
