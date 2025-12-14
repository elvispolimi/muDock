#pragma once

#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  struct queue_xsimd: public queue_cpp {
    queue_xsimd(const int _id): queue_cpp(_id) {};
  };
} // namespace mudock
