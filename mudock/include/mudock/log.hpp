#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace mudock {

  namespace log_details {
    inline void log_add_line(std::ostringstream& s) { s << std::endl; }
    template<class T, class... Ts>
    void log_add_line(std::ostringstream& s, T&& what, Ts&&... remainder) {
      s << what;
      log_add_line(s, remainder...);
    }

    class timer {
      static std::chrono::steady_clock::time_point start;

    public:
      static inline float get() {
        return std::chrono::duration<float>(std::chrono::steady_clock::now() - start).count();
      }
    };

    template<class... Ts>
    void log(Ts&&... args) {
      // declare the string stream and line composer for our log function
      std::ostringstream stream;

      // start by printing the elapsed time
      stream << '[' << std::fixed << std::setprecision(2) << std::setw(12) << std::setfill(' ')
             << timer::get() << " ] ";
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

} // namespace mudock
