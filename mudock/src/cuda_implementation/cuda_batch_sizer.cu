#include <mudock/cuda_implementation/cuda_batch_sizer.cuh>

namespace mudock {
  int compute_batch_size([[maybe_unused]] const int num_atoms, [[maybe_unused]] const int num_rotamers) {
    // TODO make something useful
    return 10;
  }
} // namespace mudock
