#pragma once

/**
 * @file buffer_alpaka.hpp
 * @brief Memory buffer specialization for the Alpaka backend.
 * @details Implements muDock's generic `buffer_impl` template for `queue_alpaka`,
 *          managing mirrored host (STL container) and device (`alpaka::Buf`) memory
 *          buffers with synchronization and copy primitives.
 */

#include <alpaka/alpaka.hpp>

#include <cassert>

#include <cstdint>
#include <memory>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/compute/buffer.hpp>
#include <optional>
#include <vector>

namespace mudock {
  /**
   * @brief Specialization of muDock's memory buffer abstraction for Alpaka.
   * @details Manages host-allocated storage (in a user-selected container, e.g. std::vector)
   *          alongside a lazily allocated device buffer (`alpaka::Buf`). Provides utilities
   *          for bidirectional memory transfers (host <-> device, device <-> device) using Alpaka views.
   *
   * @tparam container_type Host container template (e.g. std::vector).
   * @tparam T Element type stored in the buffer.
   * @tparam args Variadic template arguments for the host container allocator/traits.
   */
  template<template<class...> class container_type, typename T, class... args>
  struct buffer_impl<container_type, T, queue_alpaka, args...> {
  private:
    using dim         = queue_alpaka::dim;
    using idx         = queue_alpaka::idx;
    using dev_acc     = queue_alpaka::dev_acc;
    using buffer_type = alpaka::Buf<dev_acc, T, dim, idx>;
    using byte_type   = std::uint8_t;

    container_type<T, args...> host;              ///< Host-side memory storage container.
    std::shared_ptr<queue_alpaka> q;             ///< Shared pointer to the queue managing operations.
    std::optional<buffer_type> buffer;           ///< Optional holding the device-allocated Alpaka buffer.
    T* ptr                 = nullptr;            ///< Raw pointer to the device memory buffer.
    std::size_t size       = 0;                  ///< Current active element count.
    std::size_t alloc_size = 0;                  ///< Capacity of device allocation in elements.
    bool valid             = false;              ///< Cache coherence/validity state flag.

    /**
     * @brief Creates a 1D Alpaka extent vector for a given element count.
     * @param num_elements Number of elements.
     * @return 1D Alpaka extent vector.
     */
    static auto extent(const std::size_t num_elements) {
      return alpaka::Vec<dim, idx>{static_cast<idx>(num_elements)};
    }

    /**
     * @brief Obtains the CPU host device object used for memory views.
     * @return Host CPU Alpaka device.
     */
    static auto host_device() { return alpaka::getDevByIdx(alpaka::PlatformCpu{}, 0); }

    /**
     * @brief Allocates or expands device memory if required capacity exceeds current size.
     * @param num_elements Required capacity in elements.
     */
    void alloc_device(const std::size_t num_elements) {
      if (num_elements > alloc_size) {
        q->synchronize();
        buffer.emplace(alpaka::allocBuf<T, idx>(q->native_device(), extent(num_elements)));
        ptr        = alpaka::getPtrNative(*buffer);
        alloc_size = num_elements;
      }
      size = num_elements;
    }

  public:
    /**
     * @brief Constructs an Alpaka buffer wrapper.
     * @param _queue Shared pointer to the Alpaka queue.
     * @param num_elements Initial element capacity to allocate (default 0).
     */
    buffer_impl(std::shared_ptr<queue_alpaka> _queue, const std::size_t num_elements = 0): q(_queue) {
      if (num_elements)
        this->alloc(num_elements);
    }
    buffer_impl(const buffer_impl& other) = delete;

    buffer_impl& operator=(const buffer_impl& other) = delete;

    /// @brief Move constructor.
    buffer_impl(buffer_impl&&) = default;

    /// @brief Move assignment operator.
    buffer_impl& operator=(buffer_impl&&) = default;

    /// @brief Checks whether the device buffer currently holds valid data.
    bool is_valid() { return valid; }

    /// @brief Marks the device buffer as containing valid data.
    void set_valid() { valid = true; }

    /// @brief Marks the device buffer as invalid/stale.
    void set_not_valid() { valid = false; }

    /// @brief Destructor. Synchronizes the queue before destroying the device buffer.
    ~buffer_impl() {
      if (buffer) {
        q->synchronize();
      }
    }

    /**
     * @brief Returns a const pointer to the host memory container.
     * @return Const raw pointer to host data.
     */
    [[nodiscard]] inline auto host_pointer() const { return host.data(); }

    /**
     * @brief Returns a mutable pointer to the host memory container.
     * @return Mutable raw pointer to host data.
     */
    [[nodiscard]] inline auto host_pointer() { return host.data(); }

    /**
     * @brief Queries the number of elements allocated in the host container.
     * @return Number of elements.
     */
    [[nodiscard]] inline auto num_elements() const { return host.size(); };

    /**
     * @brief Function call operator returning mutable pointer to host data.
     * @return Mutable raw pointer to host data.
     */
    auto* operator()() { return host.data(); }

    /**
     * @brief Resizes host and device buffers to accommodate the requested number of elements.
     * @param num_elements Number of elements to allocate.
     */
    inline void alloc(const std::size_t num_elements) {
      assert(num_elements > 0 && num_elements < (std::size_t(1) << 32));
      host.resize(num_elements);
      alloc_device(num_elements);
    };

    /**
     * @brief Resizes and initializes host and device buffers to a specified byte value.
     * @param num_elements Number of elements to allocate.
     * @param value Byte value used to fill the buffers.
     */
    inline void alloc(const std::size_t num_elements, const int value) {
      assert(num_elements > 0 && num_elements < (std::size_t(1) << 32));
      host.resize(num_elements, value);
      alloc_device(num_elements);

      const auto bytes = size * sizeof(T);
      auto bytes_host  = std::vector<byte_type>(bytes, static_cast<byte_type>(value));
      auto host_view   = alpaka::createView(host_device(),
                                          bytes_host.data(),
                                          alpaka::Vec<dim, idx>{static_cast<idx>(bytes)});
      auto device_view = alpaka::createView(q->native_device(),
                                            reinterpret_cast<byte_type*>(ptr),
                                            alpaka::Vec<dim, idx>{static_cast<idx>(bytes)});
      alpaka::memcpy(q->native_queue(), device_view, host_view);
      q->synchronize();
    };

    /**
     * @brief Asynchronously copies data from host container to device buffer.
     * @param copy_size Optional number of elements to copy; if 0, copies entire host buffer.
     */
    inline void copy_host2device(const std::size_t copy_size = 0) {
      valid        = true;
      const auto n = copy_size ? copy_size : host.size();
      alloc_device(n);
      auto host_view   = alpaka::createView(host_device(), host.data(), extent(n));
      auto device_view = alpaka::createView(q->native_device(), ptr, extent(n));
      alpaka::memcpy(q->native_queue(), device_view, host_view);
    };

    /// @brief Asynchronously copies data from device buffer back to host container.
    inline void copy_device2host() {
      valid            = false;
      auto device_view = alpaka::createView(q->native_device(), ptr, extent(size));
      auto host_view   = alpaka::createView(host_device(), host.data(), extent(size));
      alpaka::memcpy(q->native_queue(), host_view, device_view);
    };

    /**
     * @brief Copies data directly from another Alpaka device buffer to this one.
     * @param other Source buffer to copy from.
     * @param n Number of elements to copy (0 = all elements in source).
     */
    inline void copy_device2device(const buffer_impl<container_type, T, queue_alpaka, args...>& other,
                                   const int n = 0) {
      valid = true;
      alloc_device(other.num_elements());
      host.resize(other.num_elements());

      const auto copy_size = n ? static_cast<std::size_t>(n) : size;
      auto src_view        = alpaka::createView(q->native_device(), other.ptr, extent(copy_size));
      auto dest_view       = alpaka::createView(q->native_device(), ptr, extent(copy_size));
      alpaka::memcpy(q->native_queue(), dest_view, src_view);
    }

    /**
     * @brief Accesses the raw device memory pointer.
     * @return Raw pointer to memory on the accelerator device.
     */
    [[nodiscard]] inline auto dev_pointer() { return ptr; }

    /**
     * @brief Returns a pointer to the internal raw device memory pointer.
     * @return Pointer to device memory pointer.
     */
    [[nodiscard]] T** dev_pointer_ref() { return &ptr; }
  };
} // namespace mudock
