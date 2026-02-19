#pragma once

#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <vector>

namespace mudock {

  template<class T>
  class safe_stack {
  public:
    using value_type = T;

  private:
    std::vector<std::unique_ptr<value_type>> stack;
    mutable std::mutex mutex;
    std::condition_variable cv;
    bool closed = false;

  public:
    // non-blocking dequeue, return nullptr if the stack is empty
    [[nodiscard]] inline auto dequeue() {
      auto new_element = std::unique_ptr<value_type>{};
      {
        std::lock_guard lock{mutex};
        if (!stack.empty()) {
          new_element = std::move(stack.back());
          stack.pop_back();
        }
      }
      return new_element;
    }

    // blocking dequeue, return nullptr if the stack is empty and closed
    std::unique_ptr<value_type> dequeue_wait() {
      std::unique_lock lock{mutex};
      // release the lock and reacquire that based on the condition
      cv.wait(lock, [&]{ return closed || !stack.empty(); });
      if (stack.empty()) return {};
      auto p = std::move(stack.back());
      stack.pop_back();
      return p;
    }

    void enqueue(std::unique_ptr<value_type> p) {
      if (!p) return;
      {
        std::lock_guard lock{mutex};
        if (closed) return;
        stack.emplace_back(std::move(p));
      }
      cv.notify_one();
    }

    // close the stack, i.e., no more element will be enqueued
    void close() {
      {
        std::lock_guard lock{mutex};
        closed = true;
      }
      cv.notify_all();
    }

    inline std::size_t clear() {
      std::lock_guard lock{mutex};
      const std::size_t n = stack.size();
      stack.clear();
      return n;
    }

    [[nodiscard]] inline std::size_t size() const {
      std::lock_guard lock{mutex};
      return stack.size();
    }
  };

} // namespace mudock
