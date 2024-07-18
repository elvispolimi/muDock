#pragma once

#include <functional>
#include <memory>
#include <thread>
#include <vector>

namespace mudock {

  class worker_interface {
  public:
    virtual void main() = 0;

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
    ~threadpool() {
      for (auto& thread: threads) {
        if (thread.joinable()) {
          thread.join();
        }
      }
    }

    template<class worker_type, class... T>
    inline void add_worker(T&&... args) {
      workers.emplace_back(std::make_unique<worker_type>(args...));
      auto* pointer = workers.back().get();
      threads.emplace_back(std::thread([pointer]() { pointer->main(); }));
    }
  };

} // namespace mudock
