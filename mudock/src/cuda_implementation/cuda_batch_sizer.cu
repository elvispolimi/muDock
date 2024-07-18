#include <mudock/cuda_implementation/cuda_batch_sizer.cuh>

namespace mudock {
  std::size_t compute_batch_size([[maybe_unused]] const std::size_t num_atoms,
                                 [[maybe_unused]] const std::size_t num_rotamers) {
    // TODO make something useful
    return 10;
  }
} // namespace mudock
