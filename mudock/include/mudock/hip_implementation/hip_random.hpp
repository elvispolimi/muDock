#pragma once

#include <hiprand/hiprand.h>
#include <hiprand/hiprand_kernel.h>
#include <mudock/hip_implementation/hip_object.hpp>

namespace mudock {
  struct hip_random_object: private hip_object<hiprandState> {
    hip_random_object()                                     = default;
    hip_random_object(const hip_random_object &)            = delete;
    hip_random_object(hip_random_object &&)                 = default;
    hip_random_object &operator=(const hip_random_object &) = delete;
    hip_random_object &operator=(hip_random_object &&)      = default;

    void alloc(const std::size_t num_elements);
    [[nodiscard]] inline auto dev_pointer() const { return hip_object<hiprandState>::dev_pointer(); }
  };

} // namespace mudock
