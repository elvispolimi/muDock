#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/log.hpp>
#include <mudock/sycl_implementation/queue_sycl.hpp>
#include <mudock/sycl_implementation/queue_sycl_impl.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <sycl/sycl.hpp>

namespace mudock {

  queue_sycl::queue_sycl(const int _id, const device_type dev_type)
      : queue(_id, dev_type), impl_(std::make_unique<impl>(_id, dev_type)) {
    assert((dev_type == device_type::GPU || dev_type == device_type::CPU) &&
           "SYCL supports only CPUs or GPUs devices");
  };
  queue_sycl::~queue_sycl() = default; // unique_ptr will destroy Impl

  queue_sycl::queue_sycl(queue_sycl&&) noexcept = default;

  void queue_sycl::alloc(void** ptr, size_t bytes) {
    queue_sycl::free(ptr);
    *ptr = sycl::malloc_device(bytes, impl_->get_queue());
    if (!*ptr)
      throw std::bad_alloc{};
  }

  void queue_sycl::free(void** ptr) {
    if (*ptr != nullptr) {
      sycl::free(*ptr, impl_->get_queue());
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
  int queue_sycl::get_preferred_workgroup_size() { return impl_->get_preferred_wg_size(); }
} // namespace mudock
