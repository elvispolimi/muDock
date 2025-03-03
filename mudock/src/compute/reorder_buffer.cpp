#include <algorithm>
#include <mudock/compute/reorder_buffer.hpp>

namespace mudock {
  reorder_buffer::reorder_buffer(std::function<int(const int, const int)> get_size) {
    // make sure to have a sizer function
    if (!get_size) {
      get_size = [](const int, const int) -> int { return int{1}; };
    }

    // populate the array of maximum sizes
    auto index = std::size_t{0};
    for (const auto& max_atom_value: atoms_clusters) {
      for (const auto& max_rotamers_value: rotamer_clusters) {
        // populate the max sizes (using the function)
        max_sizes[index] = std::min(get_size(max_atom_value, max_rotamers_value), batch::max_batch_size);

        // populate the clusters' size
        auto& cluster              = clusters[index]; // get a ref
        cluster.batch_max_atoms    = max_atom_value;
        cluster.batch_max_rotamers = max_rotamers_value;

        // update the global index on the clusters
        ++index;
      }
    }
  }

  std::pair<batch, bool> reorder_buffer::flush_one() {
    std::lock_guard lock{mutex};
    for (auto& cluster: clusters) {
      const std::size_t num_ligands = cluster.num_ligands;
      if (num_ligands > std::size_t{0}) {
        return std::make_pair(std::move(cluster), true);
      }
    }
    return std::make_pair(batch{}, false);
  }
} // namespace mudock
