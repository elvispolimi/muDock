#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/hip_implementation/geom_transform_hip.hpp>
#include <mudock/hip_implementation/mutate.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  template<>
  void geometric<queue_hip>::operator()() {
    void* args[] = {(void*) &this->args.chromosomes_per_ligand,
                    (void*) &this->args.batch_atom_stride,
                    (void*) &this->args.original_x,
                    (void*) &this->args.original_y,
                    (void*) &this->args.original_z,
                    (void*) &this->args.x_scratch,
                    (void*) &this->args.y_scratch,
                    (void*) &this->args.z_scratch,
                    (void*) &this->args.chromosomes,
                    (void*) &this->args.fragments,
                    (void*) &this->args.ligand_fragments_start,
                    (void*) &this->args.fragments_start_index,
                    (void*) &this->args.fragments_stop_index,
                    (void*) &this->args.frag_indices_start,
                    (void*) &this->args.num_rotamers,
                    (void*) &this->args.num_atoms};
    //TODO check grid/dimensions
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          this->get_queue()->launch_kernel((void*) apply_hip_kernel<max_atoms>,
                                           args,
                                           this->args.batch_ligands,
                                           BLOCK_SIZE);
        },
        this->args.batch_atom_stride,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

} // namespace mudock
