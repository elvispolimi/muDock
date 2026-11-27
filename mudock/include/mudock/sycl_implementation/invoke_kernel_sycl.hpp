#pragma once

#include <mudock/sycl_implementation/sycl_kernel_lock.hpp>
#include <mudock/sycl_implementation/queue_sycl_impl.hpp>
#include <cassert>
#include <mutex>
#include <sycl/sycl.hpp>
#include <tuple>
#include <type_traits>

namespace mudock {
  template<class F, class... Args>
  inline void queue_sycl::invoke_kernel(index3D gridDim, index3D blockDim, Args&&... args) {
    static_assert(std::is_invocable_r_v<void, F, sycl::nd_item<3>, std::decay_t<Args>...>,
                  "Kernel::operator() must be callable as: void operator()(sycl::nd_item<3>, Args...)");

    assert(gridDim.size_x() > 0 && blockDim.size_x() > 0);
    assert(gridDim.size_y() > 0 && blockDim.size_y() > 0);
    assert(gridDim.size_z() > 0 && blockDim.size_z() > 0);

    sycl::range<3> local{1, 1, static_cast<size_t>(blockDim.size_x())};
    sycl::range<3> groups{1, 1, static_cast<size_t>(gridDim.size_x())};
    sycl::range<3> global = groups * local;
    sycl::nd_range<3> nd{global, local};

    // sycl::queue q{pick_device(0, device_type::GPU)};
    sycl::event evt{};
#ifdef MUDOCK_KERNEL_LOCK
    auto* lock = get_sycl_kernel_lock(this->id);
    std::unique_lock<std::mutex> guard(lock->mutex);
    const bool has_previous_event = lock->has_event;
    const sycl::event previous_event = lock->event;
    evt = impl_->get_queue().submit([&](sycl::handler& h) {
      if (has_previous_event) {
        h.depends_on(previous_event);
      }
#else
    evt = impl_->get_queue().submit([&](sycl::handler& h) {
#endif
      F kernel{};
      const auto args_copy = std::tuple<std::decay_t<Args>...>{static_cast<std::decay_t<Args>>(args)...};
      h.parallel_for<F>(nd, [=](sycl::nd_item<3> it) {
        std::apply([&](auto... a) { kernel(it, a...); }, args_copy);
      });
    });
#ifdef MUDOCK_KERNEL_LOCK
    lock->event     = evt;
    lock->has_event = true;
#endif
  }

  template<class F, class... Args>
  inline void queue_sycl::invoke_kernel(const int gridDim, const int blockDim, Args&&... args) {
    queue_sycl::invoke_kernel<F>(index3D{gridDim, 1, 1}, index3D{blockDim, 1, 1}, args...);
  }

} // namespace mudock
