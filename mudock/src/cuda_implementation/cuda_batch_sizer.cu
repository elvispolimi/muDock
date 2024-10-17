#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cuda_implementation/cuda_batch_sizer.cuh>
#include <mudock/cuda_implementation/evaluate_fitness.cuh>
#include <mudock/utils.hpp>
#include <stdexcept>
#include <stdio.h>

#define BUCKET_MULTIPLIER 3

namespace mudock {
  int compute_batch_size(const int num_atoms, const int num_rotamers) {
    // populate the bucket dimension
    int bucket_size{0};
    constexpr_for<0, reorder_buffer::atoms_clusters.size(), 1>([&](const auto atoms_index) {
      constexpr_for<0, reorder_buffer::rotamer_clusters.size(), 1>([&](const auto rotamers_index) {
        if (num_atoms == reorder_buffer::atoms_clusters[atoms_index] &&
            num_rotamers == reorder_buffer::rotamer_clusters[rotamers_index])
          bucket_size = get_evaluate_fitness_batch<reorder_buffer::atoms_clusters[atoms_index],
                                                   reorder_buffer::rotamer_clusters[rotamers_index]>();
      });
    });
    printf("%d %d %d\n", num_atoms, num_rotamers, bucket_size);
    // TODO check if it can be made a compile error
    if (bucket_size == 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms and rotamers number which it is not handled.");
    return bucket_size * BUCKET_MULTIPLIER;
  }
} // namespace mudock
