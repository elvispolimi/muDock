#include <cassert>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/hip_implementation/geom_transform_hip.hpp>
#include <mudock/hip_implementation/hip_utils.hpp>
#include <mudock/log.hpp>
#include <mudock/hip_implementation/mutate.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  template<int MAX_ATOMS>
  batch_multiple get_geom_apply_batch(const int device_id) {
    return get_kernel_batch_multiple_hip<apply_hip<MAX_ATOMS>>(device_id,
                                                               BLOCK_SIZE,
                                                               0,
                                                               "geometric::apply_hip");
  }

  template<>
  batch_multiple get_geom_transform_batch_multiple<queue_hip>(const int atoms, std::shared_ptr<queue_hip> q_b) {
    batch_multiple bucket_multiple{};
    const int device_id = q_b->get_id();
    constexpr_for<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>([&](const auto atom_index) {
      const auto n_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
      if (atoms == n_atoms)
        bucket_multiple = get_geom_apply_batch<n_atoms>(device_id);
    });
    if (bucket_multiple.total_multiple() <= 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");
    return normalize_batch_multiple(bucket_multiple);
  }

  template<>
  void geom_kernel<queue_hip>::operator()() {
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
          q->launch_kernel((void*) apply_hip<max_atoms>, args, batch_ligands, BLOCK_SIZE);
        },
        batch_atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

} // namespace mudock
