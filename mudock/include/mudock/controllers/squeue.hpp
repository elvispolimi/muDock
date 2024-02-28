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
  // Each application's worker dock and score in batches of get_max_work() ligands
  constexpr std::size_t get_max_work() { return std::size_t{1000}; };

  // Is a Safe queue implementation
  // Enqueueing and dequeueing elements is thread safe
  template<typename T>
  class squeue {
  public:
    squeue()  = default;
    ~squeue() = default;

    squeue(const squeue&)            = delete;
    squeue(squeue&&)                 = delete;
    squeue& operator=(const squeue&) = delete;
    squeue& operator=(squeue&&)      = delete;

    // It returns ONE element of the queue if available
    // SIDE-EFFECT: the function give you ownership of the returned unique pointer
    [[nodiscard]] std::unique_ptr<T> dequeue() {
      std::lock_guard<std::mutex> guard(this->m_t);
      return dequeue_i();
    };

    // It returns X elements of the queue, where X is the minimum between the size requested and the elements still in the queue
    // SIDE-EFFECT: the function give you ownership of the returned unique pointers
    [[nodiscard]] std::vector<std::unique_ptr<T>> dequeue(const std::size_t size) {
      std::lock_guard<std::mutex> guard(this->m_t);
      std::vector<std::unique_ptr<T>> work;
      for (std::size_t index = 0; index < std::min(size, q_data.size()); index++) {
        work.emplace_back(dequeue_i());
      }
      return work;
    };

    // It enqueues ONE element in the queue
    void enqueue(std::unique_ptr<T> data) {
      std::lock_guard<std::mutex> guard(this->m_t);
      enqueue_i(std::move(data));
    };

    // It enqueues all elements in the span
    // SIDE-EFFECT: all unique ptr in the span will be nullptr upon returning
    //              the queue take ownership of the unique pointer in the span
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

    [[nodiscard]] inline std::unique_ptr<T> dequeue_i() {
      std::unique_ptr<T> res{nullptr};
      if (!q_data.empty()) {
        res = std::move(q_data.front());
        q_data.pop();
      }
      return res;
    }
  };
} // namespace mudock
