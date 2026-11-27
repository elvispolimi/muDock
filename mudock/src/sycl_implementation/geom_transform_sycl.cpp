#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/sycl_implementation/geom_transform_sycl.hpp>
#include <mudock/sycl_implementation/invoke_kernel_sycl.hpp>
#include <mudock/sycl_implementation/mutate.hpp>
#include <mudock/utils.hpp>

#ifndef MUDOCK_SYCL_WG_SIZE
  #define MUDOCK_SYCL_WG_SIZE 32
#endif

namespace mudock {
  template<>
  void geometric<queue_sycl>::operator()() {
    //TODO check grid/dimensions
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          this->get_queue()->invoke_kernel<apply_sycl<max_atoms>>(this->args.batch_ligands,
                                                  MUDOCK_SYCL_WG_SIZE,
                                                  this->args.chromosomes_per_ligand,
                                                  this->args.batch_atom_stride,
                                                  this->args.original_x,
                                                  this->args.original_y,
                                                  this->args.original_z,
                                                  this->args.x_scratch,
                                                  this->args.y_scratch,
                                                  this->args.z_scratch,
                                                  this->args.chromosomes,
                                                  this->args.fragments,
                                                  this->args.ligand_fragments_start,
                                                  this->args.fragments_start_index,
                                                  this->args.fragments_stop_index,
                                                  this->args.frag_indices_start,
                                                  this->args.num_rotamers,
                                                  this->args.num_atoms);
        },
        this->args.batch_atom_stride,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

} // namespace mudock
