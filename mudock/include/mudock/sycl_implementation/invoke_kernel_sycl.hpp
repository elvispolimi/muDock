#pragma once

#include <mudock/sycl_implementation/sycl_kernel_lock.hpp>
#include <mudock/sycl_implementation/queue_sycl_impl.hpp>
#include <algorithm>
#include <limits>
#include <mutex>
#include <string>
#include <sycl/sycl.hpp>
#include <stdexcept>
#include <typeinfo>
#include <unordered_map>

namespace mudock {
  template<class F, class... Args>
  inline void queue_sycl::invoke_kernel(index3D gridDim, index3D blockDim, Args&&... args) {
    static_assert(std::is_invocable_r_v<void, F, sycl::nd_item<3>, std::decay_t<Args>...>,
                  "Kernel::operator() must be callable as: void operator()(sycl::nd_item<3>, Args...)");

    assert(gridDim.size_x() > 0 && blockDim.size_x() > 0);
    assert(gridDim.size_y() > 0 && blockDim.size_y() > 0);
    assert(gridDim.size_z() > 0 && blockDim.size_z() > 0);

    sycl::range<3> local{(size_t) blockDim.size_x(), (size_t) blockDim.size_y(), (size_t) blockDim.size_z()};
    sycl::range<3> groups{(size_t) gridDim.size_x(), (size_t) gridDim.size_y(), (size_t) gridDim.size_z()};
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

  template<class kernel_name>
  batch_multiple queue_sycl::get_batch_multiple() {
    const sycl::queue queue = impl_->get_queue();
    const sycl::device dev = queue.get_device();
    const sycl::context ctx = queue.get_context();
    const auto dev_name = dev.get_info<sycl::info::device::name>();
    const auto dev_vendor = dev.get_info<sycl::info::device::vendor>();
    const auto dev_driver = dev.get_info<sycl::info::device::driver_version>();
    const auto kernel_key = std::string(typeid(kernel_name).name());
    const std::size_t wg_size = std::max<std::size_t>(1, MUDOCK_SYCL_WG_SIZE);
    const std::size_t dynamic_local_memory_size = 0;
    const std::string cache_key = dev_name + "|" + dev_vendor + "|" + dev_driver +
                                  "|" + kernel_key + "|" + std::to_string(wg_size);

    static std::mutex cache_mutex;
    static std::unordered_map<std::string, batch_multiple> cache_by_key;
    {
      const std::lock_guard<std::mutex> lock(cache_mutex);
      const auto it = cache_by_key.find(cache_key);
      if (it != cache_by_key.end()) {
        return it->second;
      }
    }

    const auto kid         = sycl::get_kernel_id<kernel_name>();
    const auto kb          = sycl::get_kernel_bundle<sycl::bundle_state::executable>(ctx, {dev}, {kid});
    const sycl::kernel krn = kb.get_kernel(kid);

    const std::size_t kernel_max_wg_size =
        std::max<std::size_t>(1, krn.get_info<sycl::info::kernel_device_specific::work_group_size>(dev));
    if (wg_size > kernel_max_wg_size) {
      throw std::runtime_error(
          "MUDOCK_SYCL_WG_SIZE=" + std::to_string(wg_size) +
          " exceeds the kernel maximum work-group size=" +
          std::to_string(kernel_max_wg_size));
    }

    // This is an experimental oneAPI extension. It is intentionally not
    // replaced by kernel_max_wg_size, which is only a work-group-size limit
    // and does not describe resident work-group occupancy.
    namespace syclex = sycl::ext::oneapi::experimental;
    using max_num_work_groups =
        syclex::info::kernel_queue_specific::max_num_work_groups;
    const std::size_t total_active_work_groups =
        krn.ext_oneapi_get_info<max_num_work_groups>(
            queue,
            sycl::range<1>{wg_size},
            dynamic_local_memory_size);
    if (total_active_work_groups == 0) {
      throw std::runtime_error(
          "The SYCL kernel has no resident work-group for work-group size=" +
          std::to_string(wg_size));
    }
    if (total_active_work_groups >
        static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error(
          "SYCL total active work-groups do not fit in int");
    }

    // SYCL reports a device-wide total; do not reinterpret it as
    // active work-groups per compute unit.
    const batch_multiple value{static_cast<int>(total_active_work_groups), 1};
    {
      const std::lock_guard<std::mutex> lock(cache_mutex);
      cache_by_key.emplace(cache_key, value);
    }
    return value;
  }
} // namespace mudock
