#pragma once

#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  struct queue_xsimd: public queue_cpp {
    queue_xsimd(const int _id, const device_type d_t): queue_cpp(_id, d_t) {};
  };
} // namespace mudock
