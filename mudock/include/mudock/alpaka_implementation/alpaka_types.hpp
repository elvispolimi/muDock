#pragma once

/**
 * @file alpaka_types.hpp
 * @brief Core type definitions and accelerator selections for the Alpaka backend.
 * @details Maps compile-time target backend macros (Serial, OpenMP, CUDA, HIP, SYCL)
 *          to concrete Alpaka accelerator, device, queue, and event types. Also provides
 *          primitives for synchronization and kernel locking.
 */

#include <alpaka/alpaka.hpp>
#include <mutex>
#include <memory>
#include <cstddef>
#include <type_traits>

namespace mudock::alpaka_backend {
  /// @brief One-dimensional execution dimensionality used across muDock Alpaka kernels.
  using dim = alpaka::DimInt<1u>;

  /// @brief Standard index and size type for Alpaka dimensions and buffers.
  using idx = std::size_t;

#if defined(MUDOCK_ALPAKA_BACKEND_SERIAL)
  /// @brief Accelerator type: Serial execution on single CPU thread.
  using acc = alpaka::AccCpuSerial<dim, idx>;
#elif defined(MUDOCK_ALPAKA_BACKEND_THREADS)
  /// @brief Accelerator type: CPU execution using standard C++ threads.
  using acc = alpaka::AccCpuThreads<dim, idx>;
#elif defined(MUDOCK_ALPAKA_BACKEND_TBB)
  /// @brief Accelerator type: CPU execution parallelized with Intel TBB blocks.
  using acc = alpaka::AccCpuTbbBlocks<dim, idx>;
#elif defined(MUDOCK_ALPAKA_BACKEND_OMP2)
  /// @brief Accelerator type: CPU execution parallelized with OpenMP 2.0 blocks.
  using acc = alpaka::AccCpuOmp2Blocks<dim, idx>;
#elif defined(MUDOCK_ALPAKA_BACKEND_CUDA)
  /// @brief Accelerator type: GPU execution targeting NVIDIA CUDA runtime.
  using acc = alpaka::AccGpuCudaRt<dim, idx>;
#elif defined(MUDOCK_ALPAKA_BACKEND_HIP)
  /// @brief Accelerator type: GPU execution targeting AMD HIP runtime.
  using acc = alpaka::AccGpuHipRt<dim, idx>;
#elif defined(MUDOCK_ALPAKA_BACKEND_SYCL)
  /// @brief Accelerator type: GPU execution targeting Intel SYCL runtime.
  using acc = alpaka::AccGpuSyclIntel<dim, idx>;
#else
  #error "MUDOCK_ALPAKA_BACKEND_* compile definition is required when MUDOCK_USE_ALPAKA is enabled"
#endif

  /// @brief Device type associated with the active accelerator.
  using dev_acc    = alpaka::Dev<acc>;

  /**
   * @brief Execution queue mode.
   * @details Blocking queue for CPU devices to avoid busy-waiting, NonBlocking for GPU accelerators.
   */
  using queue_mode = std::conditional_t<std::is_same_v<dev_acc, alpaka::DevCpu>,
                                        alpaka::Blocking,
                                        alpaka::NonBlocking>;

  /// @brief Command execution queue type for the configured accelerator and mode.
  using queue_acc  = alpaka::Queue<acc, queue_mode>;

  /// @brief Event type used for stream synchronization on the active queue.
  using event_acc  = alpaka::Event<queue_acc>;

  /**
   * @brief Synchronization structure for device kernel locking.
   * @details Used when concurrent kernel submission needs to be serialized across
   *          multiple host worker threads sharing the same physical device.
   */
  struct device_kernel_lock {
    std::mutex mutex;                     ///< Mutex protecting kernel submission and event recording.
    std::unique_ptr<event_acc> event;      ///< Alpaka event recorded after kernel enqueue.
    bool has_previous_event{false};       ///< Flag indicating whether a prior event is available to wait on.
  };

  /**
   * @brief Retrieves the singleton kernel lock instance for a given device.
   * @param dev_id Integer identifier of the accelerator device.
   * @param dev Reference to the Alpaka device object.
   * @return Pointer to the device_kernel_lock instance managing the specified device.
   */
  device_kernel_lock* get_kernel_lock(int dev_id, const dev_acc& dev);
} // namespace mudock::alpaka_backend

/**
 * @brief Checks if kernel locking is enabled at compile time.
 * @return True if MUDOCK_KERNEL_LOCK macro is defined, false otherwise.
 */
constexpr bool is_kernel_lock_enabled() {
#ifdef MUDOCK_KERNEL_LOCK
  return true;
#else
  return false;
#endif
}
