#include <algorithm>
#include <mudock/parallel/reorder_buffer.hpp>

namespace mudock {
  reorder_buffer::reorder_buffer(std::function<std::size_t(const std::size_t, const std::size_t)> get_size) {
    // make sure to have a sizer function
    if (!get_size) {
      get_size = [](const std::size_t, const std::size_t) -> std::size_t { return std::size_t{1}; };
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

  void reorder_buffer::flush(std::vector<batch> input_batches) {
    std::lock_guard lock{mutex};
    for (auto& cluster: clusters) {
      const auto num_ligands = cluster.num_ligands;
      if (num_ligands > std::size_t{0}) {
        auto& new_batch       = input_batches.emplace_back(batch{});
        new_batch.num_ligands = num_ligands;
        for (std::size_t i{0}; i < num_ligands; ++i) {
          new_batch.molecules[i] = std::move(cluster.molecules[i]);
        }
        cluster.num_ligands = std::size_t{0};
      }
    }
  }
} // namespace mudock
