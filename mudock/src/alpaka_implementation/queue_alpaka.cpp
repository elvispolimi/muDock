/**
 * @file queue_alpaka.cpp
 * @brief Implementation of the Alpaka compute queue and device synchronization.
 * @details Contains PIMPL definition for `queue_alpaka::impl`, device acquisition,
 *          kernel lock management, and stream synchronization wrappers.
 */

#include <alpaka/alpaka.hpp>
#include <cassert>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/log.hpp>
#include <mutex>
#include <stdexcept>

namespace mudock {
  namespace alpaka_backend {
    namespace {
      constexpr int k_max_devices = 16;
      device_kernel_lock alpaka_kernel_locks[k_max_devices];
    } // namespace

    /**
     * @brief Retrieves the device-specific kernel lock structure.
     * @details Maps CPU devices to index 0 (as only a single DevCpu device exists)
     *          and GPU devices to their corresponding index modulo `k_max_devices`.
     *          Lazily instantiates the associated Alpaka event on first access.
     * @param dev_id Requested device index.
     * @param dev Reference to the underlying Alpaka device object.
     * @return Pointer to the static `device_kernel_lock` entry.
     */
    device_kernel_lock* get_kernel_lock(int dev_id, const dev_acc& dev) {
      const int safe_id = std::is_same_v<dev_acc, alpaka::DevCpu> ? 0 : (dev_id % k_max_devices);
      if (!alpaka_kernel_locks[safe_id].event) {
        alpaka_kernel_locks[safe_id].event = std::make_unique<event_acc>(dev);
      }
      return &alpaka_kernel_locks[safe_id];
    }
  } // namespace alpaka_backend

  /**
   * @struct queue_alpaka::impl
   * @brief Private implementation holding Alpaka device and queue instances.
   * @details Isolates the heavy Alpaka template types from external translation units.
   */
  struct queue_alpaka::impl {
    dev_acc device;  ///< Native Alpaka device object.
    queue_acc queue; ///< Native Alpaka command execution queue.

    /**
     * @brief Constructs implementation by resolving the device and initializing its queue.
     * @param id Device index (automatically clamped to 0 for CPU).
     */
    impl(const int id)
        : device(
              alpaka::getDevByIdx(alpaka::Platform<acc>{}, std::is_same_v<dev_acc, alpaka::DevCpu> ? 0 : id)),
          queue(device) {}
  };

  queue_alpaka::queue_alpaka(const int _id, const device_type _dev_type)
      : queue(_id, _dev_type), impl_(std::make_unique<impl>(_id)) {
    assert((dev_type == device_type::CPU || dev_type == device_type::GPU) &&
           "Alpaka supports only CPU or GPU device tags");
  }

  queue_alpaka::~queue_alpaka() = default;

  queue_alpaka::queue_alpaka(queue_alpaka&&) noexcept = default;

  queue_alpaka::queue_acc& queue_alpaka::native_queue() { return impl_->queue; }

  const queue_alpaka::dev_acc& queue_alpaka::native_device() const { return impl_->device; }

  void queue_alpaka::alloc(void** ptr, const size_t bytes) {
    (void) ptr;
    (void) bytes;
    throw std::runtime_error("queue_alpaka::alloc is handled by buffer_impl<..., queue_alpaka>");
  }

  void queue_alpaka::free(void** ptr) {
    (void) ptr;
    throw std::runtime_error("queue_alpaka::free is handled by buffer_impl<..., queue_alpaka>");
  }

  void queue_alpaka::set_to_value(void* ptr, const size_t num_bytes, const char value) {
    (void) ptr;
    (void) num_bytes;
    (void) value;
    throw std::runtime_error("queue_alpaka::set_to_value is handled by buffer_impl<..., queue_alpaka>");
  }

  void queue_alpaka::copy_host2device(const void* host, void* device, const size_t num_bytes) {
    (void) host;
    (void) device;
    (void) num_bytes;
    throw std::runtime_error("queue_alpaka::copy_host2device is handled by buffer_impl<..., queue_alpaka>");
  }

  void queue_alpaka::copy_device2host(const void* device, void* host, const size_t num_bytes) {
    (void) device;
    (void) host;
    (void) num_bytes;
    throw std::runtime_error("queue_alpaka::copy_device2host is handled by buffer_impl<..., queue_alpaka>");
  }

  void queue_alpaka::copy_device2device(const void* device_src, void* device_dest, const size_t num_bytes) {
    (void) device_src;
    (void) device_dest;
    (void) num_bytes;
    throw std::runtime_error("queue_alpaka::copy_device2device is handled by buffer_impl<..., queue_alpaka>");
  }

  void queue_alpaka::operator()() {
    mudock::info("Worker ALPAKA on duty for device ", id, " using ", alpaka::getAccName<acc>());
  }

  void queue_alpaka::synchronize() { alpaka::wait(impl_->queue); }
} // namespace mudock
