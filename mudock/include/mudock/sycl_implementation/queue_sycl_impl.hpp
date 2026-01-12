#pragma once
#include <mudock/sycl_implementation/queue_sycl.hpp>
#include <mudock/sycl_implementation/sycl_utils.hpp>
#include <sycl/sycl.hpp>

namespace mudock {
  struct queue_sycl::impl {
    explicit impl(const int device_id, const device_type dev_type)
        : q(pick_device(device_id, dev_type), sycl::property::queue::in_order{}), d(q.get_device()) {
      preferred_wg_size = d.template get_info<sycl::info::device::sub_group_sizes>().at(0);
    };
    // no copying
    impl(const impl&)            = delete;
    impl& operator=(const impl&) = delete;

    // no moving
    impl(impl&&)            = delete;
    impl& operator=(impl&&) = delete;

    ~impl() = default;

    sycl::device get_device() { return d; }
    sycl::queue get_queue() { return q; };
    int get_preferred_wg_size() { return preferred_wg_size; }

  private:
    sycl::queue q;
    sycl::device d;
    int preferred_wg_size;
  };

} // namespace mudock
