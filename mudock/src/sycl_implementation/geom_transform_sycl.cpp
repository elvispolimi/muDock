#include <cassert>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/sycl_implementation/geom_transform_sycl.hpp>
#include <mudock/sycl_implementation/invoke_kernel_sycl.hpp>
#include <mudock/sycl_implementation/mutate.hpp>
#include <mudock/sycl_implementation/sycl_utils.hpp>
#include <mudock/log.hpp>
#include <mudock/utils.hpp>
#include <stdexcept>

#ifndef MUDOCK_SYCL_WG_SIZE
  #define MUDOCK_SYCL_WG_SIZE 32
#endif

namespace mudock {
  template<>
  batch_multiple get_geom_transform_batch_multiple<queue_sycl>(const int atoms, std::shared_ptr<queue_sycl> q_b) {
    batch_multiple bucket_multiple{};
    constexpr_for<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>([&](const auto atom_index) {
      const auto n_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
      if (atoms == n_atoms)
        bucket_multiple = get_kernel_batch_multiple_sycl<apply_sycl<n_atoms>>(q_b, "geometric::apply_sycl");
    });
    if (bucket_multiple.total_multiple() <= 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");
    return normalize_batch_multiple(bucket_multiple);
  }

  template<>
  void geom_kernel<queue_sycl>::operator()() {
    //TODO check grid/dimensions
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->invoke_kernel<apply_sycl<max_atoms>>(batch_ligands,
                                                  MUDOCK_SYCL_WG_SIZE,
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
