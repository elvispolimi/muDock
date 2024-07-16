#pragma once

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
    std::mutex mutex;

  public:
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

    inline void enqueue(std::unique_ptr<value_type> new_element) {
      if (new_element) {
        std::lock_guard lock{mutex};
        stack.emplace_back(std::move(new_element));
      }
    }
  };

} // namespace mudock
