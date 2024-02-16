#pragma once

#include <array>
#include <memory>
#include <mudock/controllers/squeue.hpp>
#include <mudock/controllers/worker.hpp>
#include <mudock/devices.hpp>
#include <mudock/molecule.hpp>
#include <span>
#include <thread>
#include <vector>

namespace mudock {
  template<language_types l_t>
  class manager {
  public:
    manager(std::shared_ptr<squeue<static_molecule>>,
            std::shared_ptr<squeue<static_molecule>>,
            const device_conf&);
    ~manager() { this->has_finished(); };

    manager(const manager&)            = delete;
    manager(manager&&)                 = delete;
    manager& operator=(const manager&) = delete;
    manager& operator=(manager&&)      = delete;

    inline void has_finished() {
      for (std::thread& w: this->workers) w.join();
      return;
    };

  protected:
    std::vector<std::thread> workers;
  };
} // namespace mudock
