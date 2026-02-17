#pragma once

#include <mutex>
#include <sycl/sycl.hpp>

namespace mudock {
#ifdef MUDOCK_KERNEL_LOCK
  struct sycl_device_kernel_lock {
    std::mutex mutex;
    sycl::event event;
    bool has_event{false};
  };

  sycl_device_kernel_lock* get_sycl_kernel_lock(int dev);
#endif
} // namespace mudock
