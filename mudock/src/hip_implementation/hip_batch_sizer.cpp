#include <mudock/hip_implementation/hip_batch_sizer.hpp>
#include <hip/hip_runtime.h>
#include <mudock/hip_implementation/evaluate_fitness.hpp>
#include <mudock/hip_implementation/hip_utils.hpp>
#include <mudock/compute/reorder_buffer.hpp>

#if defined(__HIP_PLATFORM_NVCC__) || defined(__NVCC__)
  #define BLOCK_SIZE 32
#elif defined(__HIP_PLATFORM_AMD__)
  #define BLOCK_SIZE 64
#endif

#define BUCKET_MULTIPLIER 3

namespace mudock {
  template<int MAX_ATOMS>
  int get_evaluate_fitness_batch() {
    int device_id = 0;
    MUDOCK_CHECK(hipGetDevice(&device_id));

    hipDeviceProp_t props;
    MUDOCK_CHECK(hipGetDeviceProperties(&props, device_id));

    int num_block_per_SM = 0;
    MUDOCK_CHECK(
        hipOccupancyMaxActiveBlocksPerMultiprocessor(&num_block_per_SM,
                                                     reinterpret_cast<void*>(evaluate_fitness<MAX_ATOMS>),
                                                     BLOCK_SIZE,
                                                     0));

    return num_block_per_SM * props.multiProcessorCount;
  }

  int compute_batch_size(const int num_atoms) {
    int bucket_size{0};

    constexpr_for<0, reorder_buffer::atoms_clusters.size(), 1>([&](const auto atoms_index) {
      const auto n_atoms = reorder_buffer::atoms_clusters[atoms_index];
      if (num_atoms == n_atoms)
        bucket_size = get_evaluate_fitness_batch<n_atoms>();
    });

    if (bucket_size == 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms and rotamers number which it is not handled.");

    return bucket_size * BUCKET_MULTIPLIER;
  }
} // namespace mudock
