#include <cassert>
#include <mudock/compute/queue.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/hip_implementation/adt_score_hip.hpp>
#include <mudock/hip_implementation/queue_hip.hpp>
#include <mudock/molecule.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  template<int MAX_ATOMS>
  __global__ void calc_energy(const adt_score_args args) {
    const int ligand_id = blockIdx.x;
    const int atom_stride = args.batch_atom_stride;
    const int scores_per_ligand = args.scores_per_ligand;
    fp_type* const scores_l = args.scores + ligand_id * scores_per_ligand;

    assert(blockDim.x == BLOCK_SIZE && warpSize == BLOCK_SIZE &&
           "Warpsize and the number of thread per block does not coincide");

    for (int scores_index = 0; scores_index < scores_per_ligand; ++scores_index) {
      const int coordinate_offset = (ligand_id * scores_per_ligand + scores_index) * atom_stride;
      auto pose_args = args;
      pose_args.ligand_index = ligand_id;
      pose_args.x = args.x + coordinate_offset;
      pose_args.y = args.y + coordinate_offset;
      pose_args.z = args.z + coordinate_offset;

      const fp_type score = adt_score_device_kernel::template score<MAX_ATOMS>(pose_args);
      if (threadIdx.x == 0)
        scores_l[scores_index] = score;
    }
  }

  template<>
  void adt_score<queue_hip>::operator()() {
    const auto& stage_args = this->args;
    auto q = this->get_queue();
    auto kernel_args = stage_args;
    void* args[] = {(void*) &kernel_args};
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          constexpr auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->launch_kernel((void*) calc_energy<max_atoms>, args, kernel_args.batch_ligands, BLOCK_SIZE);
        },
        kernel_args.batch_atom_stride,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

} // namespace mudock
