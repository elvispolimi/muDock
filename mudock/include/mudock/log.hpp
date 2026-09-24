#pragma once

#include <iomanip>
#include <iostream>
#include <sstream>
#include <time.h>
#include <utility>

namespace mudock {

  namespace log_details {
    inline void log_add_line(std::ostringstream& s) { s << std::endl; }
    template<class T, class... Ts>
    void log_add_line(std::ostringstream& s, T&& what, Ts&&... remainder) {
      s << what;
      log_add_line(s, remainder...);
    }

    class timer {
      static timespec start;

    public:
      static inline float get() {
        timespec current{};
        clock_gettime(CLOCK_MONOTONIC, &current);
        const auto seconds = current.tv_sec - start.tv_sec;
        const auto nanoseconds = current.tv_nsec - start.tv_nsec;
        return static_cast<float>(seconds) + static_cast<float>(nanoseconds) * 1.0e-9F;
      }
    };

    template<class... Ts>
    void log(Ts&&... args) {
      // declare the string stream and line composer for our log function
      std::ostringstream stream;

      // Prefix each line with the elapsed wall-clock time since logger startup.
      stream << '[' << std::fixed << std::setprecision(2) << std::setw(12)
             << std::setfill(' ') << timer::get() << " ] ";
      log_add_line(stream, args...);

      // print in output the line
      // NOTE: since it's a single line, the log should be thread safe
      // NOTE2: we the stderr since stdout is for results
      std::cerr << stream.str() << std::flush;
    }
  } // namespace log_details

  template<class... Ts>
  void info(Ts&&... args) {
    log_details::log("   INFO ", args...);
  }

  template<class... Ts>
  void warning(Ts&&... args) {
    log_details::log("WARNING ", args...);
  }

  template<class... Ts>
  void error(Ts&&... args) {
    log_details::log("  ERROR ", args...);
  }

  template<class... Ts>
  void stage_bucket_trace(Ts&&... args) {
#ifdef MUDOCK_ENABLE_STAGE_BUCKET_TRACE
    info(std::forward<Ts>(args)...);
#else
    (void) sizeof...(args);
#endif
  }

} // namespace mudock
