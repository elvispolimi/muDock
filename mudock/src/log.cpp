#include <chrono>
#include <mudock/log.hpp>

namespace mudock {

  namespace log_details {

    std::chrono::steady_clock::time_point timer::start = std::chrono::steady_clock::now();

  }
} // namespace mudock
