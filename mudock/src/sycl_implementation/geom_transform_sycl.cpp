#include <cassert>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/sycl_implementation/geom_transform_sycl.hpp>
#include <mudock/sycl_implementation/invoke_kernel_sycl.hpp>
#include <mudock/sycl_implementation/mutate.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  template<>
  void geom_kernel<queue_sycl>::operator()() {
    //TODO check grid/dimensions
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->invoke_kernel<apply_sycl<max_atoms>>(batch_ligands,
                                                  q->get_preferred_workgroup_size(),
                                                  chromsomes_per_ligand,
                                                  batch_atoms,
                                                  x_coords_b,
                                                  y_coords_b,
                                                  z_coords_b,
                                                  x_scratch_b,
                                                  y_scratch_b,
                                                  z_scratch_b,
                                                  chromosomes_b,
                                                  ligand_fragments_b,
                                                  ligand_fragments_start_b,
                                                  frag_start_indices_b,
                                                  frag_stop_indices_b,
                                                  frag_indices_start_b,
                                                  num_rotamers_b,
                                                  num_atoms_b);
        },
        batch_atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

} // namespace mudock
