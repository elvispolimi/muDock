#include <mudock/alpaka_implementation/queue_alpaka.hpp>

#include <alpaka/alpaka.hpp>

#include <cassert>
#include <stdexcept>

namespace mudock {
  struct queue_alpaka::impl {
    dev_acc device;
    queue_acc queue;

    impl(const int id)
        : device(alpaka::getDevByIdx(alpaka::Platform<acc>{}, id)),
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
    throw std::runtime_error("queue_alpaka::alloc is handled by object<T, queue_alpaka>");
  }

  void queue_alpaka::free(void** ptr) {
    (void) ptr;
    throw std::runtime_error("queue_alpaka::free is handled by object<T, queue_alpaka>");
  }

  void queue_alpaka::set_to_value(void* ptr, const size_t num_bytes, const char value) {
    (void) ptr;
    (void) num_bytes;
    (void) value;
    throw std::runtime_error("queue_alpaka::set_to_value is handled by object<T, queue_alpaka>");
  }

  void queue_alpaka::copy_host2device(const void* host, void* device, const size_t num_bytes) {
    (void) host;
    (void) device;
    (void) num_bytes;
    throw std::runtime_error("queue_alpaka::copy_host2device is handled by object<T, queue_alpaka>");
  }

  void queue_alpaka::copy_device2host(const void* device, void* host, const size_t num_bytes) {
    (void) device;
    (void) host;
    (void) num_bytes;
    throw std::runtime_error("queue_alpaka::copy_device2host is handled by object<T, queue_alpaka>");
  }

  void queue_alpaka::copy_device2device(const void* device_src, void* device_dest, const size_t num_bytes) {
    (void) device_src;
    (void) device_dest;
    (void) num_bytes;
    throw std::runtime_error("queue_alpaka::copy_device2device is handled by object<T, queue_alpaka>");
  }

  void queue_alpaka::operator()() {
    mudock::info("Worker ALPAKA on duty for device ", id, " using ", alpaka::getAccName<acc>());
  }

  void queue_alpaka::synchronize() { alpaka::wait(impl_->queue); }
} // namespace mudock
