#pragma once

#include <concepts>
#include <mudock/batch.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  template<typename queue_type, buffer_data_type bdt_x, buffer_data_type bdt_y, buffer_data_type bdt_z>
    requires std::derived_from<queue_type, queue>
  bool load_coordinates(batch<static_molecule> &batch,
                        std::shared_ptr<scratchpad<queue_type>> scratch,
                        const int scores_per_ligand = 1) {
    const auto batch_ligands      = batch.num_ligands;
    const auto batch_atoms        = batch.batch_max_atoms;
    const auto tot_atoms_in_batch = batch_ligands * batch_atoms;

    auto &scratch_x = (*scratch).template get<bdt_x>();
    auto &scratch_y = (*scratch).template get<bdt_y>();
    auto &scratch_z = (*scratch).template get<bdt_z>();

    if (!scratch_x.is_valid()) {
      const auto total =
          static_cast<std::size_t>(tot_atoms_in_batch) * static_cast<std::size_t>(scores_per_ligand);
      scratch_x.alloc(total);
      scratch_y.alloc(total);
      scratch_z.alloc(total);
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand = *batch.molecules[ligand_index];

        const int num_atoms = ligand.num_atoms();

        const auto x = ligand.x(), y = ligand.y(), z = ligand.z();
        for (int score_index = 0; score_index < scores_per_ligand; ++score_index) {
          const int score_offset = ligand_index * batch_atoms * scores_per_ligand +
                                   score_index * batch_atoms;
          std::memcpy((void *) (scratch_x() + score_offset), x, num_atoms * sizeof(fp_type));
          std::memcpy((void *) (scratch_y() + score_offset), y, num_atoms * sizeof(fp_type));
          std::memcpy((void *) (scratch_z() + score_offset), z, num_atoms * sizeof(fp_type));
        }
      }
      scratch_x.copy_host2device();
      scratch_y.copy_host2device();
      scratch_z.copy_host2device();
      return true;
    } else
      return false;
  };

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  bool load_coords(batch<static_molecule> &batch, std::shared_ptr<scratchpad<queue_type>> scratch) {
    return load_coordinates<queue_type,
                            buffer_data_type::X_COORDS,
                            buffer_data_type::Y_COORDS,
                            buffer_data_type::Z_COORDS>(batch, scratch, 1);
  }

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  bool load_scratchs(batch<static_molecule> &batch,
                     std::shared_ptr<scratchpad<queue_type>> scratch,
                     const int scores_per_ligand = 1) {
    return load_coordinates<queue_type,
                            buffer_data_type::X_SCRATCH,
                            buffer_data_type::Y_SCRATCH,
                            buffer_data_type::Z_SCRATCH>(batch, scratch, scores_per_ligand);
  }

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  bool load_num_rotamers(batch<static_molecule> &batch, std::shared_ptr<scratchpad<queue_type>> scratch) {
    const auto batch_ligands = batch.num_ligands;

    auto &num_rotamers_b = (*scratch).template get<buffer_data_type::NUM_ROTAMERS>();

    if (!num_rotamers_b.is_valid()) {
      num_rotamers_b.alloc(batch_ligands);
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand = *batch.molecules[ligand_index];

        num_rotamers_b()[ligand_index] = static_cast<int>(ligand.num_rotamers());
      }
      num_rotamers_b.copy_host2device();
      return true;
    } else
      return false;
  };

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  bool load_num_atoms(batch<static_molecule> &batch, std::shared_ptr<scratchpad<queue_type>> scratch) {
    const auto batch_ligands = batch.num_ligands;

    auto &num_atoms_b = (*scratch).template get<buffer_data_type::NUM_ATOMS>();

    if (!num_atoms_b.is_valid()) {
      num_atoms_b.alloc(batch_ligands);
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand = *batch.molecules[ligand_index];

        num_atoms_b()[ligand_index] = ligand.num_atoms();
      }
      num_atoms_b.copy_host2device();
      return true;
    } else
      return false;
  };
} // namespace mudock
