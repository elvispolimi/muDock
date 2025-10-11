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
#ifdef MUDOCK_BUCKET_LARGE
    // the description of how we generate the clusters
    static constexpr std::array<int, 6> atoms_clusters    = {{32, 64, 128, 160, 192, 256}};
    static constexpr std::array<int, 9> non_bond_clusters = {
        {1024, 2048, 4096, 6144, 8192, 10240, 12288, 16384, 20480}};
#else
    // the description of how we generate the clusters
    static constexpr std::array<int, 1> atoms_clusters    = {{256}};
    static constexpr std::array<int, 1> non_bond_clusters = {{20480}};
#endif
    static constexpr int get_num_atom_clusters() { return atoms_clusters.size(); };
    static constexpr int get_num_non_bond_clusters() { return non_bond_clusters.size(); };

  private:
    // the actual containers of ligand batches, with the related maximum sizes
    std::array<int, atoms_clusters.size()> max_sizes;
    std::array<batch<T>, atoms_clusters.size()> clusters;
    std::mutex mutex;

    // helper functor that given a random ligand, it will find the index of its cluster
    static constexpr auto get_flattened_index(const int num_atoms, const int num_non_bond) {
      auto index_atoms = static_cast<std::size_t>(
          std::count_if(std::begin(atoms_clusters), std::end(atoms_clusters), [&num_atoms](const auto a) {
            return a <= num_atoms;
          }));
      auto index_non_bond = static_cast<std::size_t>(
          std::count_if(std::begin(non_bond_clusters),
                        std::end(non_bond_clusters),
                        [&num_non_bond](const auto a) { return a <= num_non_bond; }));
      if (index_atoms >= get_num_atom_clusters()) {
        throw std::runtime_error("Molecule with " + std::to_string(num_atoms) + " atoms, it is too large");
      } else if (index_non_bond >= get_num_non_bond_clusters()) {
        throw std::runtime_error("Molecule with " + std::to_string(num_non_bond) +
                                 " non bonds, it is too large");
      }
      return index_atoms + index_non_bond * get_num_atom_clusters();
    }

  public:
    // the constructor will initialize the max_sizes array. The input is a function that given the number of
    // atoms and rotamers, will provide the batch size
    reorder_buffer(std::function<int(const int, const int)> get_size);

    // add the molecule to a batch. If the batch is full, return it
    std::pair<batch<T>, bool> add_ligand(std::unique_ptr<T> new_molecule) {
      std::lock_guard lock{mutex};
      const auto cluster_index =
          get_flattened_index(new_molecule->num_atoms(), new_molecule->non_bond_size());
      auto& cluster                          = clusters[cluster_index]; // take a ref (to update it)
      cluster.molecules[cluster.num_ligands] = std::move(new_molecule);
      ++cluster.num_ligands;
      return cluster.num_ligands < max_sizes[cluster_index] ? std::make_pair(batch<T>{}, false)
                                                            : std::make_pair(std::move(cluster), true);
    }

    // get the first half-empty butches inside this buffer
    std::pair<batch<T>, bool> flush_one();
  };

} // namespace mudock
