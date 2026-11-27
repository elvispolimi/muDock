#include <cassert>
#include <memory>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cuda_implementation/adt_score_cuda.cuh>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  template<int MAX_ATOMS>
  __global__ void calc_energy(const adt_score_args args) {
    const int ligand_id = blockIdx.x;
    const int atom_stride       = args.batch_atom_stride;
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

      const fp_type score = score_pose_cuda<MAX_ATOMS>(pose_args);
      if (threadIdx.x == 0)
        scores_l[scores_index] = score;
    }
  }

  template<>
  void adt_score<queue_cuda>::operator()() {
    const auto& stage_args = this->args;
    auto q = this->get_queue();
    void* kernel_parameters[] = {(void*) &stage_args};
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->launch_kernel((void*) calc_energy<max_atoms>, kernel_parameters, stage_args.batch_ligands, BLOCK_SIZE);
        },
        stage_args.batch_atom_stride,
        reorder_buffer<static_molecule>::atoms_clusters.data());

  }; // namespace mudock

} // namespace mudock
