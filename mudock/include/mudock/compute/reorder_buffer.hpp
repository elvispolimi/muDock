#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <mudock/log.hpp>
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
    static constexpr std::array<int, 8> rotamers_clusters = {{1, 4, 8, 16, 32, 64, 128, 256}};
    static constexpr std::size_t get_num_rotamer_clusters() { return rotamers_clusters.size(); };

  private:
    // the actual containers of ligand batches, with the related maximum sizes
    std::array<std::array<int, rotamers_clusters.size()>, atoms_clusters.size()> max_sizes;
    std::array<std::array<std::unique_ptr<batch<T>>, rotamers_clusters.size()>, atoms_clusters.size()> clusters;
    std::mutex mutex;

    static constexpr auto get_cluster_index(const int value, const auto& cluster_values, const char* label) {
      const auto it = std::find_if(
          std::begin(cluster_values), std::end(cluster_values), [&value](const auto max_value) {
            return value <= max_value;
          });
      const auto index = static_cast<std::size_t>(std::distance(std::begin(cluster_values), it));
      if (index >= cluster_values.size()) {
        throw std::runtime_error("Molecule with " + std::to_string(value) + " " + label +
                                 " exceeds the largest reorder-buffer bucket");
      }
      return index;
    }

  public:
    // The sizing callback receives the atom and rotamer bucket limits.
    reorder_buffer(std::function<int(const int, const int)> get_size = {}) {
      // make sure to have a sizer function
      if (!get_size) {
        get_size = [](const int, const int) -> int { return int{1}; };
      }

      // populate the array of maximum sizes
      for (std::size_t atom_index = 0; atom_index < atoms_clusters.size(); ++atom_index) {
        for (std::size_t rotamer_index = 0; rotamer_index < rotamers_clusters.size(); ++rotamer_index) {
          max_sizes[atom_index][rotamer_index] =
              std::min(get_size(atoms_clusters[atom_index], rotamers_clusters[rotamer_index]),
                       batch<T>::max_batch_size);
        }
      }
    }

    // add the molecule to a batch. If the batch is full, return it
    std::pair<batch<T>, bool> add_ligand(std::unique_ptr<T> new_molecule) {
      std::lock_guard lock{mutex};
      const auto atom_index = get_cluster_index(new_molecule->num_atoms(), atoms_clusters, "atoms");
      const auto rotamer_index = get_cluster_index(
          static_cast<int>(new_molecule->num_rotamers()), rotamers_clusters, "rotamers");
      auto& cluster_ptr = clusters[atom_index][rotamer_index];
      if (!cluster_ptr) {
        cluster_ptr = std::make_unique<batch<T>>();
        cluster_ptr->batch_max_atoms = atoms_clusters[atom_index];
        cluster_ptr->batch_max_rotamers = rotamers_clusters[rotamer_index];
        stage_bucket_trace("REORDER_BUFFER allocate atom_bucket<=",
                           atoms_clusters[atom_index],
                           " rotamer_bucket<=",
                           rotamers_clusters[rotamer_index],
                           " max_ligands=",
                           max_sizes[atom_index][rotamer_index],
                           " storage_slots=",
                           batch<T>::max_batch_size,
                           " batch_object_bytes=",
                           sizeof(batch<T>));
      }
      auto& cluster = *cluster_ptr;
      cluster.num_atoms += new_molecule->num_atoms();
      cluster.molecules[cluster.num_ligands] = std::move(new_molecule);
      ++cluster.num_ligands;
      if (cluster.num_ligands < max_sizes[atom_index][rotamer_index])
        return std::make_pair(batch<T>{}, false);

      stage_bucket_trace("REORDER_BUFFER emit atom_bucket<=",
                         cluster.batch_max_atoms,
                         " rotamer_bucket<=",
                         cluster.batch_max_rotamers,
                         " ligands=",
                         cluster.num_ligands,
                         " atoms=",
                         cluster.num_atoms,
                         " capacity=",
                         max_sizes[atom_index][rotamer_index]);
      return std::make_pair(std::move(cluster), true);
    }

    // get the first half-empty butches inside this buffer
    std::pair<batch<T>, bool> flush_one() {
      std::lock_guard lock{mutex};
      for (std::size_t atom_index = 0; atom_index < clusters.size(); ++atom_index) {
        for (std::size_t rotamer_index = 0; rotamer_index < clusters[atom_index].size(); ++rotamer_index) {
          auto& cluster = clusters[atom_index][rotamer_index];
          if (cluster && cluster->num_ligands > 0) {
            stage_bucket_trace("REORDER_BUFFER flush atom_bucket<=",
                               cluster->batch_max_atoms,
                               " rotamer_bucket<=",
                               cluster->batch_max_rotamers,
                               " ligands=",
                               cluster->num_ligands,
                               " atoms=",
                               cluster->num_atoms,
                               " capacity=",
                               max_sizes[atom_index][rotamer_index]);
            return std::make_pair(std::move(*cluster), true);
          }
        }
      }
      return std::make_pair(batch<T>{}, false);
    }
  };

} // namespace mudock
