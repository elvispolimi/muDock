#include <mudock/cuda_implementation/adadelta_convergence_cuda.cuh>

namespace mudock {
  __global__ void check_adadelta_convergence_gpu(const int population_number,
                                                 const int* __restrict__ converged_ligands_b,
                                                 int* __restrict__ ls_selected_b,
                                                 const fp_type* __restrict__ scores_b,
                                                 fp_type* __restrict__ best_scores_b,
                                                 int* __restrict__ best_scores_ages_b,
                                                 chromosome* __restrict__ chromosomes_b,
                                                 chromosome* __restrict__ best_chromosomes_b,
                                                 const std::size_t completed_iterations,
                                                 const bool finalize) {
    const int ligand_index = blockIdx.x;
    if (converged_ligands_b[ligand_index]) return;

    const std::size_t ligand_offset = static_cast<std::size_t>(ligand_index) * population_number;
    for (int individual_index = threadIdx.x; individual_index < population_number;
         individual_index += blockDim.x) {
      const auto index = ligand_offset + individual_index;
      if (!ls_selected_b[index]) continue;

      update_adadelta_convergence(scores_b[index],
                                  best_scores_b[index],
                                  best_scores_ages_b[index],
                                  ls_selected_b[index],
                                  chromosomes_b[index],
                                  best_chromosomes_b[index],
                                  completed_iterations,
                                  finalize);
    }
  }

  template<>
  void adadelta_convergence_kernel<queue_cuda>::operator()(std::size_t completed_iterations,
                                                          bool finalize) {
    void* args[] = {(void*) &population_number,
                    (void*) &converged_ligands_b,
                    (void*) &ls_selected_b,
                    (void*) &scores_b,
                    (void*) &best_scores_b,
                    (void*) &best_scores_ages_b,
                    (void*) &chromosomes_b,
                    (void*) &best_chromosomes_b,
                    (void*) &completed_iterations,
                    (void*) &finalize};
    q->launch_kernel((void*) check_adadelta_convergence_gpu, args, batch_ligands, BLOCK_SIZE);
  }
} // namespace mudock
