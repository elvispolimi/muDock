#pragma once

/**
 * @file queue_alpaka.hpp
 * @brief Declaration of the Alpaka-based compute queue.
 * @details Implements muDock's generic `queue` interface using the Alpaka performance portability
 *          library, providing asynchronous memory operations, kernel launching, and synchronization.
 */

#include <cstddef>
#include <memory>
#include <mudock/alpaka_implementation/alpaka_types.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/grid/mdindex.hpp>

namespace mudock {
  /**
   * @struct queue_alpaka
   * @brief Concrete implementation of the muDock compute queue interface backed by Alpaka.
   * @details Encapsulates an Alpaka device and command queue (`queue_acc`), managing memory
   *          allocations and transfers (Host-to-Device, Device-to-Host, Device-to-Device) as well as
   *          asynchronous kernel dispatches across CPU and GPU architectures.
   */
  struct queue_alpaka: queue {
    using dim       = alpaka_backend::dim;        ///< Dimensionality of execution spaces.
    using idx       = alpaka_backend::idx;        ///< Index and size type.
    using acc       = alpaka_backend::acc;        ///< Concrete Alpaka accelerator type.
    using dev_acc   = alpaka_backend::dev_acc;    ///< Associated device type.
    using queue_acc = alpaka_backend::queue_acc;  ///< Command queue type.

    /**
     * @brief Constructs an Alpaka queue.
     * @param _id Numeric device identifier.
     * @param _dev_type Architectural device category (CPU or GPU).
     */
    queue_alpaka(const int _id, const device_type _dev_type);

    /// @brief Destructor. Releases internal Alpaka queue and device contexts.
    ~queue_alpaka();

    queue_alpaka(const queue_alpaka&)            = delete;
    queue_alpaka& operator=(const queue_alpaka&) = delete;

    /// @brief Move constructor.
    queue_alpaka(queue_alpaka&&) noexcept;
    queue_alpaka& operator=(queue_alpaka&&) noexcept = delete;

    /**
     * @brief Enqueues a kernel task using a 3D grid dimension specification.
     * @tparam F Kernel function object type to invoke on the accelerator.
     * @tparam Args Argument types forwarded to the kernel invocation.
     * @param gridDim 3D grid dimensions (currently requires 1D launch).
     * @param args Arguments forwarded to the kernel function.
     */
    template<class F, class... Args>
    inline void invoke_kernel(const index3D gridDim, Args&&... args);

    /**
     * @brief Enqueues a kernel task using a scalar 1D grid dimension.
     * @tparam F Kernel function object type to invoke on the accelerator.
     * @tparam Args Argument types forwarded to the kernel invocation.
     * @param gridDim Number of blocks in the 1D grid.
     * @param args Arguments forwarded to the kernel function.
     */
    template<class F, class... Args>
    inline void invoke_kernel(const int gridDim, Args&&... args);

    /**
     * @brief Queries the block thread count configured for Alpaka kernels.
     * @return Number of threads per block (MUDOCK_ALPAKA_BLOCK_SIZE).
     */
    static constexpr int block_threads();

    /**
     * @brief Accesses the native Alpaka command queue.
     * @return Reference to the underlying Alpaka Queue.
     */
    queue_acc& native_queue();

    /**
     * @brief Accesses the native Alpaka device associated with this queue.
     * @return Const reference to the underlying Alpaka Device.
     */
    const dev_acc& native_device() const;

    /**
     * @brief Allocates raw memory on the accelerator device.
     * @param[out] ptr Pointer to store the allocated device address.
     * @param[in] bytes Allocation size in bytes.
     */
    void alloc(void** ptr, const size_t bytes) override;

    /**
     * @brief Deallocates raw memory from the accelerator device.
     * @param[in,out] ptr Pointer to the device address to free; set to nullptr upon release.
     */
    void free(void** ptr) override;

    /**
     * @brief Sets device memory buffer to a constant byte value.
     * @param[out] ptr Destination device memory buffer.
     * @param[in] bytes Size of memory to set in bytes.
     * @param[in] value Byte value to fill.
     */
    void set_to_value(void* ptr, const size_t bytes, const char value) override;

    /**
     * @brief Asynchronously copies memory from host to device.
     * @param[in] src Source pointer in host memory.
     * @param[out] dst Destination pointer in device memory.
     * @param[in] bytes Number of bytes to transfer.
     */
    void copy_host2device(const void* src, void* dst, const size_t bytes) override;

    /**
     * @brief Asynchronously copies memory from device to host.
     * @param[in] src Source pointer in device memory.
     * @param[out] dst Destination pointer in host memory.
     * @param[in] bytes Number of bytes to transfer.
     */
    void copy_device2host(const void* src, void* dst, const size_t bytes) override;

    /**
     * @brief Asynchronously copies memory between two device buffers.
     * @param[in] src Source pointer in device memory.
     * @param[out] dst Destination pointer in device memory.
     * @param[in] bytes Number of bytes to transfer.
     */
    void copy_device2device(const void* src, void* dst, const size_t bytes) override;

    /**
     * @brief Indicates if queue objects are required for this backend.
     * @return Always true for the Alpaka implementation.
     */
    bool obj_required() override { return true; }

    /**
     * @brief Determines whether stage bucket scheduling policy is respected.
     * @return True if operating on a GPU device, false for CPU.
     */
    bool honors_stage_bucket_policy() const override { return dev_type == device_type::GPU; }

    /// @brief Executes queue-level scheduling or progress callbacks.
    void operator()() override;

    /// @brief Blocks the host calling thread until all tasks enqueued in this queue have completed.
    void synchronize() override;

  private:
    struct impl;                          ///< Opaque forward declaration of PIMPL structure.
    std::unique_ptr<impl> impl_;          ///< Pointer to internal implementation storing Alpaka objects.
  };
} // namespace mudock
