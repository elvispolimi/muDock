#pragma once

/**
 * @file alpaka_random.hpp
 * @brief Pseudo-random number generator (PRNG) abstractions for Alpaka.
 * @details Implements architecture-aware PRNG selection (Philox for GPU, Mersenne Twister for CPU)
 *          and provides state allocation and initialization wrappers.
 */

#include <alpaka/alpaka.hpp>
#include <alpaka/rand/RandPhilox.hpp>
#include <alpaka/rand/RandStdLib.hpp>
#include <alpaka/platform/Traits.hpp>
#include <type_traits>
#include <cstddef>
#include <iterator>
#include <memory>
#include <mudock/alpaka_implementation/buffer_alpaka.hpp>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/alpaka_implementation/alpaka_types.hpp>

namespace mudock {

  /**
   * @struct rng_selector
   * @brief Selects the optimal pseudo-random number generator engine for an accelerator.
   * @details Default specialization chooses `alpaka::rand::Philox4x32x10` for counter-based,
   *          stateless generation on GPU platforms.
   *
   * @tparam TAcc Alpaka accelerator type.
   * @tparam Enable SFINAE enablement parameter.
   */
  template<typename TAcc, typename Enable = void>
  struct rng_selector {
      using type = alpaka::rand::Philox4x32x10;   ///< Engine state type (Philox4x32x10 for GPU).
      using local_ref = type;                    ///< Local reference or copy type.
  };

  /**
   * @brief CPU specialization of rng_selector selecting the Mersenne Twister engine.
   * @tparam TAcc Alpaka CPU accelerator type.
   */
  template<typename TAcc>
  struct rng_selector<TAcc, std::enable_if_t<std::is_same_v<alpaka::Platform<TAcc>, alpaka::PlatformCpu>>> {
      using type = alpaka::rand::engine::cpu::MersenneTwister; ///< Engine state type for CPU.
      using local_ref = type&;                                ///< Reference to local engine state.
  };

  /// @brief Active PRNG state type bound to the configured Alpaka backend accelerator.
  using alpaka_rand_state = typename rng_selector<alpaka_backend::acc>::type;

  /// @brief Local reference type for manipulating PRNG states within kernel functions.
  using alpaka_rand_local = typename rng_selector<alpaka_backend::acc>::local_ref;

  /**
   * @struct alpaka_random_object
   * @brief Container managing device-side PRNG state buffers for genetic search exploration.
   * @details Allocates and seeds per-thread or per-block random generator states on the accelerator device.
   */
  struct alpaka_random_object {
  public:
    /**
     * @brief Constructs an Alpaka random state manager bound to a queue.
     * @param q_ Shared pointer to the queue managing execution.
     */
    alpaka_random_object(std::shared_ptr<queue_alpaka> q_): q(q_), state(q) {};
    alpaka_random_object(const alpaka_random_object &)            = delete;
    alpaka_random_object(alpaka_random_object &&)                 = delete;
    alpaka_random_object &operator=(const alpaka_random_object &) = delete;
    alpaka_random_object &operator=(alpaka_random_object &&)      = delete;

    /**
     * @brief Allocates device PRNG states and seeds them with an explicit seed value.
     * @param num_elements Number of random generator states to allocate.
     * @param seed Numeric seed value used to initialize states.
     */
    void alloc(const std::size_t num_elements, const std::size_t seed);

    /**
     * @brief Allocates device PRNG states and seeds them with a default random sequence.
     * @param num_elements Number of random generator states to allocate.
     */
    void alloc(const std::size_t num_elements);

    /**
     * @brief Accesses raw pointer to device-side RNG states.
     * @return Pointer to device memory holding RNG state array.
     */
    [[nodiscard]] inline auto dev_pointer() { return state.dev_pointer(); }

    /**
     * @brief Accesses pointer to the raw device pointer variable.
     * @return Pointer to device memory pointer.
     */
    [[nodiscard]] inline alpaka_rand_state **dev_pointer_ref() { return state.dev_pointer_ref(); }

  private:
    std::shared_ptr<queue_alpaka> q;                        ///< Bound queue instance.
    buffer_vector<alpaka_rand_state, queue_alpaka> state;   ///< Device buffer vector of RNG states.
  };
} // namespace mudock
