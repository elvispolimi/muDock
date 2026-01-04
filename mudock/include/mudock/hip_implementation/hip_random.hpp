#pragma once

#include <hiprand/hiprand.h>
#include <hiprand/hiprand_kernel.h>
#include <mudock/compute/object.hpp>
#include <mudock/hip_implementation/queue_hip.hpp>

namespace mudock {
  struct hip_random_object {
    hip_random_object(std::shared_ptr<queue_hip> q_): q(q_), state(q) {};
    hip_random_object(const hip_random_object &)            = delete;
    hip_random_object(hip_random_object &&)                 = delete;
    hip_random_object &operator=(const hip_random_object &) = delete;
    hip_random_object &operator=(hip_random_object &&)      = delete;

    void alloc(const std::size_t num_elements, const std::size_t seed);
    void alloc(const std::size_t num_elements);

    [[nodiscard]] inline auto dev_pointer() const { return state.dev_pointer(); }
    [[nodiscard]] inline hiprandState **dev_pointer_ref() { return state.dev_pointer_ref(); }

  private:
    std::shared_ptr<queue_hip> q;
    object<hiprandState> state;
  };

} // namespace mudock
