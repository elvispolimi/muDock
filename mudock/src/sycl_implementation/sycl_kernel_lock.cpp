#include <functional>
#include <memory>
#include <mudock/compute/devices_memory.hpp>
#include <mudock/sycl_implementation/sycl_kernel_lock.hpp>

namespace mudock {
#ifdef MUDOCK_KERNEL_LOCK
  namespace {
    constexpr int k_max_devices_kernel_lock = 16;

    device_memory_array<k_max_devices_kernel_lock, sycl_device_kernel_lock>* get_kernel_lock_storage() {
      // Intentionally leaked to avoid static destruction after SYCL runtime teardown.
      static auto* storage = new device_memory_array<k_max_devices_kernel_lock, sycl_device_kernel_lock>();
      return storage;
    }
  } // namespace

  sycl_device_kernel_lock* get_sycl_kernel_lock(const int dev) {
    auto* storage = get_kernel_lock_storage();
    storage->init(dev, std::function<std::unique_ptr<sycl_device_kernel_lock>()>([]() {
      return std::make_unique<sycl_device_kernel_lock>();
    }));
    return storage->v[dev].get_data();
  }
#endif
} // namespace mudock
