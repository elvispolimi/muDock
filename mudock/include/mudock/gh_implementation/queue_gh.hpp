#pragma once

#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  struct queue_gh: public queue_cpp {
    queue_gh(const int _id, const device_type dev_type): queue_cpp(_id, dev_type) {};
  };
} // namespace mudock
