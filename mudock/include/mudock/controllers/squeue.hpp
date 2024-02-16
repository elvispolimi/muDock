#pragma once

#include <algorithm>
#include <cstddef>
#include <memory>
#include <mudock/type_alias.hpp>
#include <mutex>
#include <queue>
#include <span>
#include <utility>
#include <vector>

namespace mudock {
  constexpr index_type get_max_work() { return index_type{1000}; };

  template<typename T>
  class squeue {
  public:
    squeue()                         = default;
    ~squeue()                        = default;
    squeue(const squeue&)            = delete;
    squeue(squeue&&)                 = delete;
    squeue& operator=(const squeue&) = delete;
    squeue& operator=(squeue&&)      = delete;

    std::unique_ptr<T> dequeue() {
      std::lock_guard<std::mutex> guard(this->m_t);
      return dequeue_i();
    };

    std::vector<std::unique_ptr<T>> dequeue(const std::size_t size) {
      std::lock_guard<std::mutex> guard(this->m_t);
      std::vector<std::unique_ptr<T>> work;
      for (std::size_t index = 0; index < std::min(size, q_data.size()); index++) {
        work.emplace_back(dequeue_i());
      }
      return work;
    };

    void enqueue(std::unique_ptr<T> data) {
      std::lock_guard<std::mutex> guard(this->m_t);
      enqueue_i(std::move(data));
    };

    void enqueue(std::span<std::unique_ptr<T>> data) {
      std::lock_guard<std::mutex> guard(this->m_t);
      for (auto& t: data) enqueue_i(std::move(t));
    };

  private:
    std::queue<std::unique_ptr<T>> q_data;
    std::mutex m_t;

    template<typename C>
    inline void enqueue_i(C&& data) {
      q_data.emplace(std::move(data));
    }

    inline std::unique_ptr<T> dequeue_i() {
      std::unique_ptr<T> res{nullptr};
      if (!q_data.empty()) {
        res = std::move(q_data.front());
        q_data.pop();
      }
      return res;
    }
  };
} // namespace mudock
