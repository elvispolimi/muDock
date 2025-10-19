#include <algorithm>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/sycl_implementation/evaluate_fitness.hpp>
#include <mudock/sycl_implementation/sycl_batch_sizer.hpp>
#include <mudock/sycl_implementation/virtual_screen.hpp>
#include <mudock/utils.hpp>
#include <sycl/sycl.hpp>

namespace mudock {
  int compute_batch_size(const sycl::device& d, const int num_atoms, const int num_non_bonds) {
    return compute_batch_size_vs(d, num_atoms, num_non_bonds);
  }
} // namespace mudock
