#pragma once

#include <cstdint>

namespace mudock {
  std::size_t compute_batch_size(const std::size_t num_atoms, const std::size_t num_rotamers);
}
