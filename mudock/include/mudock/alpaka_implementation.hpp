#pragma once

#ifdef MUDOCK_USE_ALPAKA
  #include <mudock/alpaka_implementation/alpaka_implementation.hpp>
#else
  #include <mudock/log.hpp>
namespace mudock {
  inline void alpaka_backend_disabled() { warning("The Alpaka implementation is disabled"); }
} // namespace mudock
#endif
