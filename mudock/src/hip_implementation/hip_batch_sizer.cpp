#include <mudock/hip_implementation/hip_batch_sizer.hpp>

namespace mudock {
  int compute_batch_size([[maybe_unused]] const int num_atoms, [[maybe_unused]] const int num_rotamers) {
    // TODO make something useful
    return 1000;
  }
} // namespace mudock
