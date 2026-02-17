#pragma once

#include <mudock/compute/devices_memory.hpp>
#include <mudock/sycl_implementation/queue_sycl_impl.hpp>
#include <mutex>
#include <string>
#include <sycl/sycl.hpp>
#include <typeinfo>
#include <unordered_map>

namespace mudock {
#ifdef MUDOCK_KERNEL_LOCK
  namespace {
    constexpr int k_max_devices_kernel_lock = 16;

    struct device_kernel_lock {
      std::mutex mutex;
      sycl::event event;
      bool has_event{false};
    };

    device_memory_array<k_max_devices_kernel_lock, device_kernel_lock> sycl_kernel_locks;

    device_kernel_lock* get_kernel_lock(const int dev) {
      sycl_kernel_locks.init(dev, std::function<std::unique_ptr<device_kernel_lock>()>([]() {
        return std::make_unique<device_kernel_lock>();
      }));
      return sycl_kernel_locks.v[dev].get_data();
    }
  } // namespace
#endif

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
    auto* lock = get_kernel_lock(this->id);
    std::unique_lock<std::mutex> guard(lock->mutex);
    if (lock->has_event) {
      lock->event.wait_and_throw();
    }
    evt = impl_->get_queue().submit([&](sycl::handler& h) {
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
  int queue_sycl::get_batch_size() {
    const sycl::device dev      = impl_->get_device();
    const auto dev_name         = dev.get_info<sycl::info::device::name>();
    const auto dev_vendor       = dev.get_info<sycl::info::device::vendor>();
    const auto dev_driver       = dev.get_info<sycl::info::device::driver_version>();
    const auto kernel_key       = std::string(typeid(kernel_name).name());
    const std::string cache_key = dev_name + "|" + dev_vendor + "|" + dev_driver + "|" + kernel_key;

    static std::mutex cache_mutex;
    static std::unordered_map<std::string, int> cache_by_key;
    {
      const std::lock_guard<std::mutex> lock(cache_mutex);
      const auto it = cache_by_key.find(cache_key);
      if (it != cache_by_key.end()) {
        return it->second;
      }
    }

    const sycl::context ctx{dev};
    const int compute_units         = dev.get_info<sycl::info::device::max_compute_units>();
    const auto subgroups            = dev.get_info<sycl::info::device::sub_group_sizes>();
    const std::size_t subgroup_size = subgroups.empty() ? 1u : subgroups[0];

    // Fetch kernel device-specific limits for this device
    const auto kid         = sycl::get_kernel_id<kernel_name>();
    const auto kb          = sycl::get_kernel_bundle<sycl::bundle_state::executable>(ctx, {dev}, {kid});
    const sycl::kernel krn = kb.get_kernel(kid);

    // Maximum work-group size the device allows for this kernel
    const std::size_t max_wg_size = krn.get_info<sycl::info::kernel_device_specific::work_group_size>(dev);

    // Upper bound on concurrently resident work-groups per CU from wg-size alone
    const std::size_t wg_per_cu_cap =
        std::max<std::size_t>(1, max_wg_size / std::max<std::size_t>(1, subgroup_size));

    // Portable heuristic: try to keep 2–4 work-groups per CU if possible.
    // You can tune this number per backend/workload.
    const std::size_t target_wg_per_cu = std::min<std::size_t>(wg_per_cu_cap, 4);

    const int value = static_cast<int>(target_wg_per_cu * compute_units);
    {
      const std::lock_guard<std::mutex> lock(cache_mutex);
      cache_by_key.emplace(cache_key, value);
    }
    return value;
  }
} // namespace mudock
