#pragma once

#include <concepts>
#include <mudock/batch.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  template<typename queue_type, typename layout_set = all_scratch_layouts>
    requires std::derived_from<queue_type, queue>
  bool load_coordinates(batch<static_molecule> &batch,
                        std::shared_ptr<scratchpad<queue_type, layout_set>> scratch,
                        const int scores_per_ligand) {
    const auto batch_ligands      = batch.num_ligands;
    const auto batch_atoms        = batch.batch_max_atoms;
    const auto tot_atoms_in_batch = batch_ligands * batch_atoms;
    auto& coordinates = scratch->scratch_coordinates();

    if (!coordinates.is_valid()) {
      const auto total =
          static_cast<std::size_t>(tot_atoms_in_batch) * static_cast<std::size_t>(scores_per_ligand);
      coordinates.alloc(total);
      auto host_view = coordinates.host_view();
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand = *batch.molecules[ligand_index];

        const int num_atoms = ligand.num_atoms();

        const auto x = ligand.x(), y = ligand.y(), z = ligand.z();
        for (int score_index = 0; score_index < scores_per_ligand; ++score_index) {
          const int score_offset = ligand_index * batch_atoms * scores_per_ligand +
                                   score_index * batch_atoms;
          std::memcpy(host_view.x().data() + score_offset, x, num_atoms * sizeof(fp_type));
          std::memcpy(host_view.y().data() + score_offset, y, num_atoms * sizeof(fp_type));
          std::memcpy(host_view.z().data() + score_offset, z, num_atoms * sizeof(fp_type));
        }
      }
      coordinates.copy_host2device();
      return true;
    } else
      return false;
  };

  template<typename queue_type, typename layout_set = all_scratch_layouts>
    requires std::derived_from<queue_type, queue>
  bool load_coords(batch<static_molecule> &batch, std::shared_ptr<scratchpad<queue_type, layout_set>> scratch) {
    return load_coordinates(batch, scratch, 1);
  }

  template<typename queue_type, typename layout_set = all_scratch_layouts>
    requires std::derived_from<queue_type, queue>
  bool load_scratchs(batch<static_molecule> &batch,
                     std::shared_ptr<scratchpad<queue_type, layout_set>> scratch,
                     const int scores_per_ligand = 1) {
    return load_coordinates(batch, scratch, scores_per_ligand);
  }

  template<typename queue_type, typename layout_set = all_scratch_layouts>
    requires std::derived_from<queue_type, queue>
  bool load_num_rotamers(batch<static_molecule> &batch, std::shared_ptr<scratchpad<queue_type, layout_set>> scratch) {
    const auto batch_ligands = batch.num_ligands;

    auto &batch_info = scratch->batch_info();

    if (!batch_info.is_valid()) {
      batch_info.alloc(batch_ligands);
      auto view = batch_info.host_view();
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand = *batch.molecules[ligand_index];

        view.num_rotamers()[ligand_index] = static_cast<int>(ligand.num_rotamers());
        view.num_atoms()[ligand_index] = ligand.num_atoms();
      }
      batch_info.copy_host2device();
      return true;
    } else
      return false;
  };

  template<typename queue_type, typename layout_set = all_scratch_layouts>
    requires std::derived_from<queue_type, queue>
  bool load_num_atoms(batch<static_molecule> &batch, std::shared_ptr<scratchpad<queue_type, layout_set>> scratch) {
    const auto batch_ligands = batch.num_ligands;

    auto &batch_info = scratch->batch_info();

    if (!batch_info.is_valid()) {
      batch_info.alloc(batch_ligands);
      auto view = batch_info.host_view();
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand = *batch.molecules[ligand_index];

        view.num_atoms()[ligand_index] = ligand.num_atoms();
        view.num_rotamers()[ligand_index] = static_cast<int>(ligand.num_rotamers());
      }
      batch_info.copy_host2device();
      return true;
    } else
      return false;
  };
} // namespace mudock
