#pragma once

/**
 * @file invoke_kernel_alpaka.hpp
 * @brief Inline template implementations for launching kernels via queue_alpaka.
 * @details Implements work-division calculation, task creation, optional device-level
 *          kernel serialization locking, and queue enqueueing for Alpaka kernels.
 */

#include <cassert>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mutex>
#include <type_traits>
#include <utility>

namespace mudock {
  /**
   * @brief Enqueues an Alpaka kernel using 3D grid dimensions.
   * @details Constructs a 1D `alpaka::WorkDivMembers` from `gridDim.size_x()`, the static
   *          block dimension `block_threads()`, and 1 element per thread. Creates a task
   *          with `alpaka::createTaskKernel` and pushes it into `native_queue()`.
   *          If `is_kernel_lock_enabled()` evaluates to true, serializes launches across
   *          worker threads using a mutex and Alpaka events.
   *
   * @tparam F Function object (functor) type representing the Alpaka kernel.
   * @tparam Args Forwarded parameter types passed to the kernel's `operator()`.
   * @param gridDim 3D grid size vector (size_y and size_z must currently be 1).
   * @param args Arguments forwarded to the kernel invocation.
   */
  template<class F, class... Args>
  inline void queue_alpaka::invoke_kernel(const index3D gridDim, Args&&... args) {
    assert(gridDim.size_x() > 0);
    assert(gridDim.size_y() == 1 && "queue_alpaka currently supports only 1D launches");
    assert(gridDim.size_z() == 1 && "queue_alpaka currently supports only 1D launches");

    const auto grid  = alpaka::Vec<dim, idx>{static_cast<idx>(gridDim.size_x())};
    const auto block = alpaka::Vec<dim, idx>{static_cast<idx>(block_threads())};
    const auto elems = alpaka::Vec<dim, idx>{idx{1}};

    const auto work_div = alpaka::WorkDivMembers<dim, idx>{grid, block, elems};
    auto task           = alpaka::createTaskKernel<acc>(work_div,
                                              F{},
                                              static_cast<std::decay_t<Args>>(std::forward<Args>(args))...);
                                              
    if constexpr (is_kernel_lock_enabled()) {
      auto* lock = alpaka_backend::get_kernel_lock(id, native_device());
      std::unique_lock<std::mutex> guard(lock->mutex);
      if (lock->has_previous_event) {
        alpaka::wait(native_queue(), *(lock->event));
      }
      alpaka::enqueue(native_queue(), task);
      alpaka::enqueue(native_queue(), *(lock->event));
      lock->has_previous_event = true;
    } else {
      alpaka::enqueue(native_queue(), task);
    }
  }

  /**
   * @brief Enqueues an Alpaka kernel using a 1D scalar grid dimension.
   * @tparam F Function object (functor) type representing the Alpaka kernel.
   * @tparam Args Forwarded parameter types passed to the kernel's `operator()`.
   * @param gridDim Number of blocks in the 1D grid.
   * @param args Arguments forwarded to the kernel invocation.
   */
  template<class F, class... Args>
  inline void queue_alpaka::invoke_kernel(const int gridDim, Args&&... args) {
    assert(gridDim >= 0);
    invoke_kernel<F>(index3D{gridDim, 1, 1}, std::forward<Args>(args)...);
  }

  /**
   * @brief Returns the fixed block thread count configured for Alpaka kernels.
   * @return Number of threads per block (MUDOCK_ALPAKA_BLOCK_SIZE).
   */
  inline constexpr int queue_alpaka::block_threads() { return MUDOCK_ALPAKA_BLOCK_SIZE; }
} // namespace mudock
