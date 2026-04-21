#include <mudock/tbb_implementation/runtime_services.hpp>

#include <chrono>
#include <utility>

namespace mudock::detail {

  periodic_observer::~periodic_observer() {
    stop();
    join();
  }

  void periodic_observer::start(std::optional<double> period_sec, std::function<void()> callback) {
    stop();
    join();

    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop_requested_ = false;
    }

    if (!period_sec || *period_sec <= 0.0) {
      return;
    }

    thread_ = std::thread([this, period = *period_sec, callback = std::move(callback)]() mutable {
      while (true) {
        std::unique_lock<std::mutex> lock(mutex_);
        const bool stop =
            cv_.wait_for(lock, std::chrono::duration<double>(period), [&]() { return stop_requested_; });
        if (stop) {
          return;
        }
        lock.unlock();
        callback();
      }
    });
  }

  void periodic_observer::stop() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop_requested_ = true;
    }
    cv_.notify_one();
  }

  void periodic_observer::join() {
    if (thread_.joinable()) {
      thread_.join();
    }
  }

  deadline_timer::~deadline_timer() {
    cancel();
    join();
  }

  void deadline_timer::start(std::optional<double> time_limit_sec, std::function<void()> callback) {
    cancel();
    join();

    {
      std::lock_guard<std::mutex> lock(mutex_);
      cancelled_ = false;
    }

    if (!time_limit_sec || *time_limit_sec <= 0.0) {
      return;
    }

    thread_ = std::thread([this, limit = *time_limit_sec, callback = std::move(callback)]() mutable {
      std::unique_lock<std::mutex> lock(mutex_);
      const bool cancelled =
          cv_.wait_for(lock, std::chrono::duration<double>(limit), [&]() { return cancelled_; });
      if (cancelled) {
        return;
      }
      lock.unlock();
      callback();
    });
  }

  void deadline_timer::cancel() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      cancelled_ = true;
    }
    cv_.notify_one();
  }

  void deadline_timer::join() {
    if (thread_.joinable()) {
      thread_.join();
    }
  }

} // namespace mudock::detail
