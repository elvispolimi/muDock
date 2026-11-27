#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/devices.hpp>
#include <mudock/sycl_implementation/adt_score_sycl.hpp>
#include <mudock/sycl_implementation/invoke_kernel_sycl.hpp>
#include <mudock/sycl_implementation/queue_sycl.hpp>
#include <mudock/utils.hpp>

#ifndef MUDOCK_SYCL_WG_SIZE
  #define MUDOCK_SYCL_WG_SIZE 32
#endif

namespace mudock {

  template<>
  void adt_score<queue_sycl>::operator()() {
    const auto& stage_args = this->args;
    auto q                 = this->get_queue();
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->invoke_kernel<calc_energy<max_atoms>>(stage_args.batch_ligands, MUDOCK_SYCL_WG_SIZE, stage_args);
        },
        stage_args.batch_atom_stride,
        reorder_buffer<static_molecule>::atoms_clusters.data());

  }; // namespace mudock

} // namespace mudock
