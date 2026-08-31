#include <time.h>
#include <mudock/log.hpp>

namespace mudock {

  namespace log_details {

    timespec timer::start = [] {
      timespec value{};
      clock_gettime(CLOCK_MONOTONIC, &value);
      return value;
    }();

  } // namespace log_details
} // namespace mudock
