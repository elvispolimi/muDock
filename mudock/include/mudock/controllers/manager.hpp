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
  // The manager is responsible to spawn the application's threads, which execute the worker function (in worker.hpp)
  template<language_types l_t>
  class manager {
  public:
    // Manager constructor based on the choosen language spawns the threads which do dock and score stages
    manager(std::shared_ptr<squeue<static_molecule>>,
            std::shared_ptr<squeue<static_molecule>>,
            const device_conf&);
    // The manager constructor waits for all the thread to finish computation
    ~manager() { this->has_finished(); };

    // The manager cannot be moved not copied, only one for application is instantiated in the main
    manager(const manager&)            = delete;
    manager(manager&&)                 = delete;
    manager& operator=(const manager&) = delete;
    manager& operator=(manager&&)      = delete;

    // To check if all the manager's thread have completed
    inline void has_finished() {
      for (std::thread& w: this->workers) w.join();
      return;
    };

  protected:
    // Vector of all computational threads which perform dock and score
    std::vector<std::thread> workers;
  };
} // namespace mudock
