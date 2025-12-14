#pragma once

#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {
  template<typename T>
  using buffer_cpp = buffer_vector<T, queue_cpp>;
}
