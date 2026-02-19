#pragma once

#include <functional>
#include <memory>
#include <thread>
#include <vector>

namespace mudock {

  class worker_interface {
  public:
    virtual void main() = 0;

    // TODO ask Gadio
    virtual ~worker_interface() {}
  };

  class threadpool {
    std::vector<std::thread> threads;
    std::vector<std::unique_ptr<worker_interface>> workers;

  public:
    threadpool()                             = default;
    threadpool(threadpool&&)                 = default;
    threadpool(const threadpool&)            = delete;
    threadpool& operator=(threadpool&&)      = default;
    threadpool& operator=(const threadpool&) = delete;
    ~threadpool() { wait(); }

    template<class worker_type>
    inline void add_worker(worker_type wt) {
      // workers.emplace_back(std::make_unique<worker_type>(args...));
      // auto* pointer = workers.back().get();
      // threads.emplace_back(std::thread([pointer]() { pointer->main(); }));
      // TODO check me if it works
      // threads.emplace_back(std::thread([&wt]() mutable { wt.main(); }));
      threads.emplace_back(&worker_type::main, std::move(wt));
    }

    inline void wait() {
      for (auto& thread: threads) {
        if (thread.joinable()) {
          thread.join();
        }
      }
    }
  };

} // namespace mudock
