#pragma once

#include <alpaka/alpaka.hpp>
#include <alpaka/rand/RandPhilox.hpp>
#include <cstddef>
#include <iterator>
#include <memory>
#include <mudock/alpaka_implementation/buffer_alpaka.hpp>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>

namespace mudock {

  using alpaka_rand_state = alpaka::rand::Philox4x32x10;

  struct alpaka_random_object {
  public:
    alpaka_random_object(std::shared_ptr<queue_alpaka> q_): q(q_), state(q) {};
    alpaka_random_object(const alpaka_random_object &)            = delete;
    alpaka_random_object(alpaka_random_object &&)                 = delete;
    alpaka_random_object &operator=(const alpaka_random_object &) = delete;
    alpaka_random_object &operator=(alpaka_random_object &&)      = delete;

    void alloc(const std::size_t num_elements, const std::size_t seed);
    void alloc(const std::size_t num_elements);

    [[nodiscard]] inline auto dev_pointer() { return state.dev_pointer(); }
    [[nodiscard]] inline alpaka_rand_state **dev_pointer_ref() { return state.dev_pointer_ref(); }

  private:
    std::shared_ptr<queue_alpaka> q;
    buffer_vector<alpaka_rand_state, queue_alpaka> state;
  };
} // namespace mudock
