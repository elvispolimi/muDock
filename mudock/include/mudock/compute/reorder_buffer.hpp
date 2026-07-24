#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <mudock/batch.hpp>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mudock {

  // this is a buffer that re-arrange the ligands to cluster them in the same batch, if they have similar
  // features, i.e. a similar number of atoms and rotamers
  template<class T>
  class reorder_buffer {
  public:
#ifdef MUDOCK_ATOM_CLUSTER_LEVEL_EXTREME
    // the description of how we generate the clusters
    static constexpr std::array<int, 9> atoms_clusters = {{16, 32, 58, 64, 96, 128, 160, 192, 256}};
#elif defined(MUDOCK_ATOM_CLUSTER_LEVEL_LARGE)
    // the description of how we generate the clusters
    static constexpr std::array<int, 6> atoms_clusters = {{32, 64, 128, 160, 192, 256}};
#elif defined(MUDOCK_ATOM_CLUSTER_LEVEL_MEDIUM)
    // the description of how we generate the clusters
    static constexpr std::array<int, 4> atoms_clusters = {{32, 64, 128, 256}};
#else
    // the description of how we generate the clusters
    static constexpr std::array<int, 1> atoms_clusters = {{256}};
#endif
    static constexpr std::size_t get_num_atom_clusters() { return atoms_clusters.size(); };

  private:
    // the actual containers of ligand batches, with the related maximum sizes
    std::array<int, atoms_clusters.size()> max_sizes;
    std::array<batch<T>, atoms_clusters.size()> clusters;
    std::mutex mutex;

    // helper functor that given a random ligand, it will find the index of its cluster
    static constexpr auto get_flattened_index(const int num_atoms) {
      const auto it = std::find_if(
          std::begin(atoms_clusters), std::end(atoms_clusters), [&num_atoms](const auto a) {
            return num_atoms <= a;
          });
      auto index_atoms = static_cast<std::size_t>(std::distance(std::begin(atoms_clusters), it));
      if (index_atoms >= get_num_atom_clusters()) {
        throw std::runtime_error("Molecule with " + std::to_string(num_atoms) + " atoms, it is too large");
      }
      return index_atoms;
    }

  public:
    // the constructor will initialize the max_sizes array. The input is a function that given the number of
    // atoms and rotamers, will provide the batch size
    reorder_buffer(std::function<int(const int)> get_size = {}) {
      // make sure to have a sizer function
      if (!get_size) {
        get_size = [](const int) -> int { return int{1}; };
      }

      // populate the array of maximum sizes
      auto index = std::size_t{0};
      for (const auto& max_atom_value: atoms_clusters) {
        // populate the max sizes (using the function)
        max_sizes[index] = std::min(get_size(max_atom_value), batch<T>::max_batch_size);

        // populate the clusters' size
        auto& cluster           = clusters[index]; // get a ref
        cluster.batch_max_atoms = max_atom_value;

        // update the global index on the clusters
        ++index;
      }
    }

    // add the molecule to a batch. If the batch is full, return it
    std::pair<batch<T>, bool> add_ligand(std::unique_ptr<T> new_molecule) {
      std::lock_guard lock{mutex};
      const auto cluster_index               = get_flattened_index(new_molecule->num_atoms());
      auto& cluster                          = clusters[cluster_index]; // take a ref (to update it)
      cluster.molecules[cluster.num_ligands] = std::move(new_molecule);
      ++cluster.num_ligands;
      return cluster.num_ligands < max_sizes[cluster_index] ? std::make_pair(batch<T>{}, false)
                                                            : std::make_pair(std::move(cluster), true);
    }

    // get the first half-empty butches inside this buffer
    std::pair<batch<T>, bool> flush_one() {
      std::lock_guard lock{mutex};
      for (auto& cluster: clusters) {
        const std::size_t num_ligands = cluster.num_ligands;
        if (num_ligands > std::size_t{0}) {
          return std::make_pair(std::move(cluster), true);
        }
      }
      return std::make_pair(batch<T>{}, false);
    }
  };

} // namespace mudock
