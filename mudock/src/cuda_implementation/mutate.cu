#include <cstddef>
#include <mudock/cuda_implementation/geometric_transformations.cuh>
#include <mudock/cuda_implementation/mutate.cuh>

namespace mudock {

  __device__ void apply_cuda(fp_type* __restrict__ x,
                             fp_type* __restrict__ y,
                             fp_type* __restrict__ z,
                             const fp_type* chromosome,
                             const int* __restrict__ fragments,
                             const std::size_t* __restrict__ fragments_start_index,
                             const std::size_t* __restrict__ fragments_stop_index,
                             const std::size_t num_rotamers,
                             const std::size_t stride_atoms,
                             const std::size_t num_atoms) {
    // apply rigid transformations
    translate_molecule(x, y, z, chromosome, chromosome + 1, chromosome + 2, num_atoms);
    rotate_molecule(x, y, z, chromosome + 3, chromosome + 4, chromosome + 5, num_atoms);

    // change the molecule shape
    for (std::size_t i = 0; i < num_rotamers; ++i) {
      const int* bitmask = fragments + i * stride_atoms;
      rotate_fragment(x,
                      y,
                      z,
                      bitmask,
                      fragments_start_index[i],
                      fragments_stop_index[i],
                      chromosome + 6 + i,
                      num_atoms);
    }
  }

} // namespace mudock
