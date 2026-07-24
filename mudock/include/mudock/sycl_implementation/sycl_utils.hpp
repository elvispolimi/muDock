#pragma once

#include <algorithm>
#include <memory>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/devices.hpp>
#include <mudock/sycl_implementation/queue_sycl.hpp>
#include <sycl/sycl.hpp>

namespace mudock {

  inline sycl::device pick_device(const int device_id, const device_type type) {
    std::vector<sycl::device> matches;

    for (const auto& dev: sycl::device::get_devices()) {
      if (type == device_type::GPU && dev.is_gpu())
        matches.push_back(dev);
      else if (type == device_type::CPU && dev.is_cpu())
        matches.push_back(dev);
    }

    assert(!matches.empty() && "SYCL no requsted device found");
    assert(device_id < static_cast<int>(matches.size()) && "SYCL requsted device ID not found");

    return matches[static_cast<size_t>(device_id)];
  }

  template<class kernel_name>
  inline batch_multiple get_kernel_batch_multiple_sycl(std::shared_ptr<queue_sycl> q_b,
                                                       const char* kernel_label = "unknown_kernel") {
    (void) kernel_label;
    return q_b->template get_batch_multiple<kernel_name>();
  }

} // namespace mudock
