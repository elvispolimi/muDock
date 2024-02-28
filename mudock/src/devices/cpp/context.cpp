#include <cstdint>
#include <mudock/devices/cpp/context.hpp>

namespace mudock {
  template<>
  device_context<language_types::CPP>::device_context(const std::size_t id) {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(id, &cpuset); // Set affinity to CPU 0
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
  };

  template<>
  device_context<language_types::CPP>::~device_context(){};

} // namespace mudock
