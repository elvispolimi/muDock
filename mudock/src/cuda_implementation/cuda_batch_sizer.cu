#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cuda_implementation/cuda_batch_sizer.cuh>
#include <mudock/cuda_implementation/evaluate_fitness.cuh>
#include <mudock/utils.hpp>
#include <stdexcept>

#define BUCKET_MULTIPLIER 3

namespace mudock {
  template<int MAX_ATOMS>
  int get_evaluate_fitness_batch() {
    int device_id = 0;
    MUDOCK_CHECK(cudaGetDevice(&device_id));
    cudaDeviceProp props;
    MUDOCK_CHECK(cudaGetDeviceProperties(&props, device_id));
    int num_block_per_SM = 0;
    MUDOCK_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&num_block_per_SM,
                                                               evaluate_fitness<MAX_ATOMS>,
                                                               BLOCK_SIZE,
                                                               0));
    // TODO check the return value
    return num_block_per_SM * props.multiProcessorCount;
  }

  int compute_batch_size(const int num_atoms) {
    // populate the bucket dimension
    int bucket_size{0};
    constexpr_for<0, reorder_buffer::atoms_clusters.size(), 1>([&](const auto atoms_index) {
      const auto n_atoms = reorder_buffer::atoms_clusters[atoms_index];
      if (num_atoms == n_atoms)
        bucket_size = get_evaluate_fitness_batch<n_atoms>();
    });
    // TODO check if it can be made a compile error
    if (bucket_size == 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms and rotamers number which it is not handled.");
    return bucket_size * BUCKET_MULTIPLIER;
  }
} // namespace mudock
