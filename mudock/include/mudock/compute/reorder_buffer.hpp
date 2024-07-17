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
  class reorder_buffer {
  public:
    // the description of how we generate the clusters
    static constexpr std::array<std::size_t, 6> atoms_clusters   = {{32, 64, 128, 160, 192, 256}};
    static constexpr std::array<std::size_t, 6> rotamer_clusters = {{1, 2, 4, 8, 16, 32}};

  private:
    // the actual containers of ligand batches, with the related maximum sizes
    std::array<batch, atoms_clusters.size() * rotamer_clusters.size()> clusters;
    std::array<std::size_t, atoms_clusters.size() * rotamer_clusters.size()> max_sizes;
    std::mutex mutex;

    // helper functor that given a random ligand, it will find the index of its cluster
    static constexpr auto get_flattened_index(const std::size_t num_atoms, const std::size_t num_rotamers) {
      auto index_atoms = static_cast<std::size_t>(
          std::count_if(std::begin(atoms_clusters), std::end(atoms_clusters), [&num_atoms](const auto a) {
            return a >= num_atoms;
          }));
      if (index_atoms >= atoms_clusters.size()) {
        throw std::runtime_error("Molecule with " + std::to_string(num_atoms) + " atoms, it is too large");
      }
      auto index_rotamers = static_cast<std::size_t>(
          std::count_if(std::begin(rotamer_clusters),
                        std::end(rotamer_clusters),
                        [&num_rotamers](const auto r) { return r >= num_rotamers; }));
      if (index_rotamers >= rotamer_clusters.size()) {
        throw std::runtime_error("Molecule with " + std::to_string(num_rotamers) +
                                 " rotamers, it is too flexible");
      }
      return (index_atoms * atoms_clusters.size()) + index_rotamers;
    }

  public:
    // the constructor will initialize the max_sizes array. The input is a function that given the number of
    // atoms and rotamers, will provide the batch size
    reorder_buffer(std::function<std::size_t(const std::size_t, const std::size_t)> get_size);

    // add the molecule to a batch. If the batch is full, return it
    template<class container_type>
      requires is_container_specification<container_type>
    std::pair<batch&&, bool> add_ligand(std::unique_ptr<molecule<container_type>> new_molecule) {
      std::lock_guard lock{mutex};
      const auto cluster_index = get_flattened_index(new_molecule->num_atoms(), new_molecule->num_rotamers());
      auto& cluster            = clusters[cluster_index]; // take a ref (to update it)
      cluster.molecules[cluster.num_ligands] = std::move(new_molecule);
      ++cluster.num_ligands;
      return cluster.num_ligands < max_sizes[cluster_index] ? std::make_pair(batch{}, false)
                                                            : std::make_pair(std::move(cluster), true);
    }

    // get the first half-empty butches inside this buffer
    std::pair<batch&&, bool> flush_one();
  };

} // namespace mudock
