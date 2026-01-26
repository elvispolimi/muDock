#pragma once
#include <mudock/sycl_implementation/queue_sycl.hpp>
#include <mudock/sycl_implementation/sycl_utils.hpp>
#include <sycl/sycl.hpp>

namespace mudock {
  struct queue_sycl::impl {
    explicit impl(const int device_id, const device_type dev_type)
        : q(pick_device(device_id, dev_type), sycl::property::queue::in_order{}), d(q.get_device()) {}
    // no copying
    impl(const impl&)            = delete;
    impl& operator=(const impl&) = delete;

    // no moving
    impl(impl&&)            = delete;
    impl& operator=(impl&&) = delete;

    ~impl() = default;

    sycl::device get_device() { return d; }
    sycl::queue get_queue() { return q; };

  private:
    sycl::queue q;
    sycl::device d;
  };

} // namespace mudock
