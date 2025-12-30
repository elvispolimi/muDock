#include <cassert>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cuda_implementation/geom_transform_cuda.cuh>
#include <mudock/cuda_implementation/mutate.cuh>
#include <mudock/utils.hpp>

namespace mudock {
  template<>
  void geom_kernel<queue_cuda>::operator()() {
    void* args[] = {(void*) &chromsomes_per_ligand,
                    (void*) &batch_atoms,
                    (void*) &x_coords_b,
                    (void*) &y_coords_b,
                    (void*) &z_coords_b,
                    (void*) &x_scratch_b,
                    (void*) &y_scratch_b,
                    (void*) &z_scratch_b,
                    (void*) &chromosomes_b,
                    (void*) &ligand_fragments_b,
                    (void*) &ligand_fragments_start_b,
                    (void*) &frag_start_indices_b,
                    (void*) &frag_stop_indices_b,
                    (void*) &frag_indices_start_b,
                    (void*) &num_rotamers_b,
                    (void*) &num_atoms_b};
    //TODO check grid/dimensions
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->launch_kernel((void*) apply_cuda<max_atoms>, args, batch_ligands);
        },
        batch_atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

} // namespace mudock
