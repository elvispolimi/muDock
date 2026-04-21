#pragma once

#include <cassert>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <memory>
#include <mudock/log.hpp>
#include <mutex>

namespace mudock {

  template<typename T>
  class safe_queue {
    // this is the actual container of the queue
    std::deque<std::unique_ptr<T>> buffer;
    std::size_t max_buffer_size;

    // this is a counter only for statistical purposes
    std::size_t global_counter;

    // tto signal when to exit
    bool signal_terminate;

    // these variables avoids busy waiting and actually prevent the thread to consume CPU cycles
    mutable std::mutex queue_mutex;
    std::condition_variable inwork_available;
    std::condition_variable outwork_available;

  public:
    using value_type     = T;
    using value_ptr_type = std::unique_ptr<T>;

    safe_queue(void): max_buffer_size(1), global_counter(0), signal_terminate(false) {}

    // this queue cannot be copied or moved around
    safe_queue(const safe_queue &) = delete;
    safe_queue(safe_queue &&)      = delete;

    inline void initialize(const std::size_t max_queue_size) {
      std::unique_lock<std::mutex> lock(queue_mutex);
      max_buffer_size = max_queue_size;
    }

    inline std::size_t size(void) const {
      std::unique_lock<std::mutex> lock(queue_mutex);
      return buffer.size();
    }

    inline std::size_t max_size(void) const { return max_buffer_size; }

    inline bool empty(void) const {
      std::unique_lock<std::mutex> lock(queue_mutex);
      return buffer.empty();
    }

    inline std::size_t get_global_counter(void) const {
      std::unique_lock<std::mutex> lock(queue_mutex);
      return global_counter;
    }

    inline void send_terminate_signal(void) {
      std::unique_lock<std::mutex> lock(queue_mutex);
      signal_terminate = true;
      inwork_available.notify_all();
      outwork_available.notify_all();
    }

    inline void clear_terminate_signal(void) {
      std::unique_lock<std::mutex> lock(queue_mutex);
      signal_terminate = false;
    }

    inline std::size_t clear() {
      std::unique_lock<std::mutex> lock(queue_mutex);
      const std::size_t n = buffer.size();
      buffer.clear();
      inwork_available.notify_all();
      return n;
    }

    inline bool get_terminate_signal(void) {
      std::unique_lock<std::mutex> lock(queue_mutex);
      return signal_terminate;
    }

    // this method attempt to remove and return an element from the queue.
    // NOTE: this method might block the execution of the application.
    value_ptr_type dequeue() {
      // wait until there is some event
      std::unique_lock<std::mutex> lock(queue_mutex);
      while (!signal_terminate && buffer.empty()) { // spourious events might happens!
        outwork_available.wait(lock);               // release lock -> wait for a wake_up -> reaquire lock
      }

      // get a new job from the input queue (if any)
      if (!buffer.empty()) {
        auto output_data = std::move(buffer.back());
        buffer.pop_back();
        inwork_available.notify_one();
        return output_data;
      } else {
        assert(signal_terminate);
        return value_ptr_type{};
      }
    }

    // this method attempt to insert an element in the queue. In case of success the queue owns the element and the user is not
    // allowed to dereference the pointer. If the operation fails, i.e. the termination signal is set, the owner
    // of the data is still the caller
    bool enqueue(value_ptr_type &input_data) {
      // try to enqueue the element (if there is enough space)
      std::unique_lock<std::mutex> lock(queue_mutex);
      while (!signal_terminate && (buffer.size() >= max_buffer_size)) { // spourious events might happens!
        inwork_available.wait(lock); // release lock -> wait for a wake_up -> reaquire lock
      }

      // if the terminate signal is set, do not enqueue the element
      if (signal_terminate) {
        return false;
      }

      // otherwise, store the data in the back of the container
      buffer.emplace_front(std::move(input_data));
      outwork_available.notify_one();
      global_counter += 1;
      return true;
    }
  };
} // namespace mudock
