#pragma once

#include <condition_variable>
#include <functional>
#include <mutex>
#include <optional>
#include <thread>

namespace mudock::detail {

  class periodic_observer {
    std::mutex mutex_;
    std::condition_variable cv_;
    bool stop_requested_ = false;
    std::thread thread_;

  public:
    periodic_observer() = default;
    ~periodic_observer();

    periodic_observer(const periodic_observer&) = delete;
    periodic_observer& operator=(const periodic_observer&) = delete;
    periodic_observer(periodic_observer&&) = delete;
    periodic_observer& operator=(periodic_observer&&) = delete;

    void start(std::optional<double> period_sec, std::function<void()> callback);
    void stop();
    void join();
  };

  class deadline_timer {
    std::mutex mutex_;
    std::condition_variable cv_;
    bool cancelled_ = false;
    std::thread thread_;

  public:
    deadline_timer() = default;
    ~deadline_timer();

    deadline_timer(const deadline_timer&) = delete;
    deadline_timer& operator=(const deadline_timer&) = delete;
    deadline_timer(deadline_timer&&) = delete;
    deadline_timer& operator=(deadline_timer&&) = delete;

    void start(std::optional<double> time_limit_sec, std::function<void()> callback);
    void cancel();
    void join();
  };

} // namespace mudock::detail
