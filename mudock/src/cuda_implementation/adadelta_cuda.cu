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
                                            int *__restrict__ active_individuals_b,    
                                            const fp_type rho,
                                            const fp_type epsilon) {

    const int ligand_id        = blockIdx.x;
    const int local_thread_id  = threadIdx.x;
    assert(blockDim.x == BLOCK_SIZE && warpSize == BLOCK_SIZE &&
           "Warpsize and the number of thread per block does not coincide");

    chromosome *__restrict__ population_l         = population_b + ligand_id * individuals_per_ligand;
    chromosome *__restrict__ e_g2_l               = adadelta_e_g2_b + ligand_id * individuals_per_ligand;
    chromosome *__restrict__ e_dw2_l              = adadelta_e_dw2_b + ligand_id * individuals_per_ligand;
    gradient   *__restrict__ gradients_l          = gradients_b + ligand_id * individuals_per_ligand;
    const int  *__restrict__ active_individuals_l = active_individuals_b + ligand_id * individuals_per_ligand;
    
    const int num_rotamers = num_rotamers_b[ligand_id];
    const int gradient_size = 6 + num_rotamers;
    
    for (int individual_index = local_thread_id; individual_index < individuals_per_ligand; individual_index += blockDim.x) {
      if (!active_individuals_l[individual_index]){
        continue;
      }
      
      gradient   &grad = gradients_l[individual_index];
      chromosome &w = population_l[individual_index];
      chromosome &E_g2_i = e_g2_l[individual_index];
      chromosome &E_dw2_i = e_dw2_l[individual_index];
      // Apply AdaDelta update for each dimension
      for (int d = 0; d < gradient_size; ++d) {
        E_g2_i[d] = rho * E_g2_i[d] + (fp_type{1} - rho) * grad[d] * grad[d];
        
        const fp_type rms_g = std::sqrt(E_g2_i[d] + epsilon);
        
        const fp_type rms_dw = std::sqrt(E_dw2_i[d] + epsilon);
        
        fp_type delta_w = -(rms_dw / rms_g) * grad[d];
        
        E_dw2_i[d] = rho * E_dw2_i[d] + (fp_type{1} - rho) * delta_w * delta_w;
        w[d] = w[d] + delta_w;
      }
      
    }

  }

  
  template<>
  void adadelta_kernel<queue_cuda>::operator()() {
    void* args[] = {(void*) &batch_ligands,
                    (void*) &individuals_per_ligand,
                    (void*) &gradients_b,
                    (void*) &population_b,
                    (void*) &num_rotamers_b,
                    (void*) &adadelta_e_g2_b,
                    (void*) &adadelta_e_dw2_b,
                    (void*) &active_individuals_b,
                    (void*) &rho,
                    (void*) &epsilon};
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
      [&](const auto atom_index) {
        const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
        q->launch_kernel((void*) apply_adadelta_update_gpu<max_atoms>, args, batch_ligands, BLOCK_SIZE);
      },
      batch_atoms,
      reorder_buffer<static_molecule>::atoms_clusters.data());

  }

} // namespace mudock