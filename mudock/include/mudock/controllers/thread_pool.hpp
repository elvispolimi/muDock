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
  // The thread pool is responsible to spawn the application's threads, which execute the worker function (in worker.hpp)
  // For ech new kernel type a new thread pool has to be provided
  template<kernel_type l_t>
  class thread_pool {
  public:
    // Manager constructor based on the choosen kernel spawns the threads which do dock and score stages
    thread_pool(std::shared_ptr<squeue<static_molecule>>,
                std::shared_ptr<squeue<static_molecule>>,
                const device_conf&);
    // The thread pool constructor waits for all the thread to finish computation
    ~thread_pool() { this->has_finished(); };

    // The thread pool cannot be moved not copied, only one for application is instantiated in the main
    thread_pool(const thread_pool&)            = delete;
    thread_pool(thread_pool&&)                 = delete;
    thread_pool& operator=(const thread_pool&) = delete;
    thread_pool& operator=(thread_pool&&)      = delete;

    // To check if all the threads in the pool have completed
    inline void has_finished() {
      for (std::thread& w: this->workers) w.join();
      return;
    };

  protected:
    // Vector of all computational threads which perform dock and score
    std::vector<std::thread> workers;
  };
} // namespace mudock
