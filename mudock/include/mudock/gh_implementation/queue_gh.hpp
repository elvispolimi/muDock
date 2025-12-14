#pragma once

#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  struct queue_gh: public queue_cpp {
    queue_gh(const int _id): queue_cpp(_id) {};
  };
} // namespace mudock
