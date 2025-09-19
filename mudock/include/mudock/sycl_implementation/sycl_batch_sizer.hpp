#pragma once
#include <cstdint>
#include <sycl/sycl.hpp>

namespace mudock {
  int compute_batch_size(const sycl::device&, const int num_atoms);
}
