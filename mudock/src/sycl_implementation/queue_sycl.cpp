#include <algorithm>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/log.hpp>
#include <mudock/sycl_implementation/queue_sycl.hpp>
#include <mudock/sycl_implementation/queue_sycl_impl.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <sycl/sycl.hpp>

namespace mudock {

  queue_sycl::queue_sycl(const int _id,
                         const device_type dev_type,
                         std::shared_ptr<device_memory_tracker> tracker)
      : queue(_id, dev_type, std::move(tracker)), impl_(std::make_unique<impl>(_id, dev_type)) {
    assert((dev_type == device_type::GPU || dev_type == device_type::CPU) &&
           "SYCL supports only CPUs or GPUs devices");
  };
  queue_sycl::~queue_sycl() = default; // unique_ptr will destroy Impl

  queue_sycl::queue_sycl(queue_sycl&&) noexcept = default;

  void queue_sycl::alloc(void** ptr, size_t bytes) {
    *ptr = sycl::malloc_device(bytes, impl_->get_queue());
    if (!*ptr)
      throw std::bad_alloc{};
    impl_->allocated_bytes += bytes;
    impl_->peak_allocated_bytes = std::max(impl_->peak_allocated_bytes, impl_->allocated_bytes);
    memory_tracker->allocate(bytes);
  }

  void queue_sycl::free(void** ptr, const size_t bytes) {
    if (*ptr != nullptr) {
      sycl::free(*ptr, impl_->get_queue());
      assert(impl_->allocated_bytes >= bytes);
      impl_->allocated_bytes -= bytes;
      memory_tracker->release(bytes);
      *ptr = nullptr;
    }
  }

  void queue_sycl::set_to_value(void* ptr, const size_t num_bytes, const char value) {
    impl_->get_queue().memset(ptr, value, num_bytes);
  }

  void queue_sycl::copy_host2device(const void* host, void* device, size_t num_bytes) {
    impl_->get_queue().memcpy(device, host, num_bytes);
  }

  void queue_sycl::copy_device2host(const void* device, void* host, size_t num_bytes) {
    impl_->get_queue().memcpy(host, device, num_bytes);
  }

  void queue_sycl::copy_device2device(const void* device_src, void* device_dest, size_t num_bytes) {
    impl_->get_queue().memcpy(device_dest, device_src, num_bytes);
  }

  void queue_sycl::synchronize() { impl_->get_queue().wait_and_throw(); }

  void queue_sycl::operator()() {
    // Setup default queue/device if you want a global singleton.
  }

  std::size_t queue_sycl::allocated_bytes() const { return impl_->allocated_bytes; }
  std::size_t queue_sycl::peak_allocated_bytes() const { return impl_->peak_allocated_bytes; }
} // namespace mudock
