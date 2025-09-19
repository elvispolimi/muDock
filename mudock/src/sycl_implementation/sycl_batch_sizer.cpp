#include <algorithm>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/sycl_implementation/evaluate_fitness.hpp>
#include <mudock/sycl_implementation/sycl_batch_sizer.hpp>
#include <mudock/utils.hpp>
#include <sycl/sycl.hpp>

#define BUCKET_MULTIPLIER 3

namespace mudock {
  template<int MAX_ATOMS>
  int get_evaluate_fitness_batch(const sycl::device& dev) {
    const sycl::context ctx{dev};
    const int compute_units = dev.get_info<sycl::info::device::max_compute_units>();
    const int subgroup_size = dev.get_info<sycl::info::device::sub_group_sizes>()[0];

    // Fetch kernel device-specific limits for this device
    const auto kid         = sycl::get_kernel_id<evaluate_fitness_kernel_tag<MAX_ATOMS>>();
    const auto kb          = sycl::get_kernel_bundle<sycl::bundle_state::executable>(ctx, {dev}, {kid});
    const sycl::kernel krn = kb.get_kernel(kid);

    // Maximum work-group size the device allows for this kernel
    const std::size_t max_wg_size = krn.get_info<sycl::info::kernel_device_specific::work_group_size>(dev);

    // Preferred multiple (warp/wavefront alignment); useful when picking BLOCK_SIZE
    const std::size_t preferred_multiple =
        krn.get_info<sycl::info::kernel_device_specific::preferred_work_group_size_multiple>(dev);

    // Upper bound on concurrently resident work-groups per CU from wg-size alone
    const std::size_t wg_per_cu_cap = std::max<std::size_t>(1, max_wg_size / subgroup_size);

    // Portable heuristic: try to keep 2–4 work-groups per CU if possible.
    // You can tune this number per backend/workload.
    const std::size_t target_wg_per_cu = std::min<std::size_t>(wg_per_cu_cap, 4);

    return static_cast<int>(target_wg_per_cu * compute_units);
  }

  int compute_batch_size(const sycl::device& d, const int num_atoms) {
    // populate the bucket dimension
    int bucket_size{0};
    constexpr_for<0, reorder_buffer::atoms_clusters.size(), 1>([&](const auto atoms_index) {
      const auto n_atoms = reorder_buffer::atoms_clusters[atoms_index];
      if (num_atoms == n_atoms)
        bucket_size = get_evaluate_fitness_batch<n_atoms>(d);
    });
    // TODO check if it can be made a compile error
    if (bucket_size == 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");
    return bucket_size * BUCKET_MULTIPLIER;
  }
} // namespace mudock
