#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <mudock/molecule.hpp>

namespace mudock {

  // this struct describes a bundle of ligands. To limit the number of dynamic mallocs, we use a static
  // allocation. For this reason we have a maximum number of ligands in a batch, and the actual size
  struct batch {
    static constexpr auto max_batch_size = std::size_t{100};

    std::array<std::unique_ptr<static_molecule>, max_batch_size> molecules;
    std::size_t num_ligands        = 0;
    std::size_t batch_max_atoms    = 0;
    std::size_t batch_max_rotamers = 0;

    batch()  = default;
    ~batch() = default;
    batch(batch&& other) {
      molecules.swap(other.molecules);
      num_ligands        = other.num_ligands;
      other.num_ligands  = 0;
      batch_max_atoms    = other.batch_max_atoms;
      batch_max_rotamers = other.batch_max_rotamers;
    }
    batch(const batch& other)      = delete;
    batch& operator=(batch&&)      = default;
    batch& operator=(const batch&) = delete;
  };

} // namespace mudock
