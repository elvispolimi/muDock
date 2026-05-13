#pragma once

#include <algorithm>
#include <cstddef>
#include <memory>
#include <mudock/batch.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<typename queue_type>
  struct stage {
    static std::size_t get_shared_ligand_mem(const int, const knobs&) { return 0; }
    static std::size_t get_private_ligand_mem(const int, const knobs&) { return 0; }
    static int get_ligand_mem(const int atoms, const knobs conf) {
      return static_cast<int>(get_shared_ligand_mem(atoms, conf) + get_private_ligand_mem(atoms, conf));
    }
    static int get_batch_multiple(const int, std::shared_ptr<queue_type>, const knobs&) { return 1; }
    static int get_batch_size(const int,
                              std::shared_ptr<queue_type>,
                              const knobs&,
                              const size_t max_bucket_size) {
      return static_cast<int>(std::max<size_t>(1, max_bucket_size));
    }

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
