#include <mudock/type_alias.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>
#include <mudock/compute/adadelta_kernel.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cuda_implementation/adadelta_cuda.cuh>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/molecule.hpp>
#include <mudock/utils.hpp>
#include <mudock/cuda_implementation/cuda_utils.cuh>
#include <cmath>

namespace mudock {
  template<int MAX_ATOMS>
  __global__ void apply_adadelta_update_gpu(const int batch_ligands,
                                            const int individuals_per_ligand,
                                            gradient *__restrict__ gradients_b,
                                            chromosome *__restrict__ population_b,
                                            int* __restrict__ num_rotamers_b,
                                            chromosome *__restrict__ adadelta_e_g2_b,
                                            chromosome *__restrict__ adadelta_e_dw2_b,
                                            int *__restrict__ ls_selected_b,
                                            const fp_type rho,
                                            const fp_type epsilon,
                                            int *__restrict__ converged_ligands_b) {

    const int ligand_id        = blockIdx.x;
    const int local_thread_id  = threadIdx.x;
    assert(blockDim.x == BLOCK_SIZE && warpSize == BLOCK_SIZE &&
           "Warpsize and the number of thread per block does not coincide");

    if (converged_ligands_b[ligand_id]) return;
    chromosome *__restrict__ population_l         = population_b + ligand_id * individuals_per_ligand;
    chromosome *__restrict__ e_g2_l               = adadelta_e_g2_b + ligand_id * individuals_per_ligand;
    chromosome *__restrict__ e_dw2_l              = adadelta_e_dw2_b + ligand_id * individuals_per_ligand;
    gradient   *__restrict__ gradients_l          = gradients_b + ligand_id * individuals_per_ligand;
    const int  *__restrict__ ls_selected_l = ls_selected_b + ligand_id * individuals_per_ligand;

    const int num_rotamers  = num_rotamers_b[ligand_id];
    const int gradient_size = 6 + num_rotamers;

    const int total = individuals_per_ligand * gradient_size;
    for (int idx = local_thread_id; idx < total; idx += blockDim.x) {
      const int individual_index = idx / gradient_size;
      const int d                = idx - individual_index * gradient_size;
      if (!ls_selected_l[individual_index]) continue;

      const fp_type g = gradients_l[individual_index][d];
      fp_type& w      = population_l[individual_index][d];
      fp_type& E_g2   = e_g2_l[individual_index][d];
      fp_type& E_dw2  = e_dw2_l[individual_index][d];

      E_g2 = rho * E_g2 + (fp_type{1} - rho) * g * g;

      const fp_type rms_g  = std::sqrt(E_g2 + epsilon);
      const fp_type rms_dw = std::sqrt(E_dw2 + epsilon);

      const fp_type delta_w = -(rms_dw / rms_g) * g;

      E_dw2 = rho * E_dw2 + (fp_type{1} - rho) * delta_w * delta_w;

      w = w + delta_w;
    }
  }

  
  template<>
  void adadelta_kernel<queue_cuda>::operator()() {
    void* args[] = {(void*) &batch_ligands,
                    (void*) &individuals_per_ligand,
                    (void*) &gradients_b,
                    (void*) &chromosomes_b,
                    (void*) &num_rotamers_b,
                    (void*) &adadelta_e_g2_b,
                    (void*) &adadelta_e_dw2_b,
                    (void*) &ls_selected_b,
                    (void*) &rho,
                    (void*) &epsilon,
                    (void*) &converged_ligands_b};
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
      [&](const auto atom_index) {
        const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
        q->launch_kernel((void*) apply_adadelta_update_gpu<max_atoms>, args, batch_ligands, BLOCK_SIZE);
      },
      batch_atoms,
      reorder_buffer<static_molecule>::atoms_clusters.data());

  }

  template<int MAX_ATOMS>
  batch_multiple get_adadelta_batch(const int device_id) {
    return get_kernel_batch_multiple_cuda<apply_adadelta_update_gpu<MAX_ATOMS>>(device_id,
                                                                  BLOCK_SIZE,
                                                                  0,
                                                                  "adadelta::apply_adadelta_update_gpu");
  }

  template<>
  batch_multiple get_adadelta_batch_multiple<queue_cuda>(const int atoms, std::shared_ptr<queue_cuda> q_b) {
    batch_multiple bucket_multiple{};
    const int device_id = q_b->get_id();
    constexpr_for<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>([&](const auto atom_index) {
      const auto n_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
      if (atoms == n_atoms)
        bucket_multiple = get_adadelta_batch<n_atoms>(device_id);
    });
    if (bucket_multiple.total_multiple() <= 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");
    return normalize_batch_multiple(bucket_multiple);
  }

} // namespace mudock