#include <mudock/cpp_implementation/adadelta_convergence_cpp.hpp>

namespace mudock {
  namespace {
    void check_adadelta_convergence(const int batch_ligands,
                                    const int individuals_per_ligand,
                                    const int* converged_ligands,
                                    int* selected,
                                    const fp_type* scores,
                                    fp_type* best_scores,
                                    int* ages,
                                    chromosome* chromosomes,
                                    chromosome* best_chromosomes,
                                    const std::size_t completed_iterations,
                                    const bool finalize) {
      for (int ligand = 0; ligand < batch_ligands; ++ligand) {
        if (converged_ligands[ligand]) continue;
        for (int individual = 0; individual < individuals_per_ligand; ++individual) {
          const auto index = static_cast<std::size_t>(ligand) * static_cast<std::size_t>(individuals_per_ligand) +
                             static_cast<std::size_t>(individual);
          update_adadelta_convergence(scores[index], best_scores[index], ages[index], selected[index],
                                      chromosomes[index], best_chromosomes[index], completed_iterations, finalize);
        }
      }
    }
  }

  template<>
  void adadelta_convergence_kernel<queue_cpp>::operator()(const std::size_t completed_iterations,
                                                         const bool finalize) {
    q->invoke_kernel<adadelta_convergence_region_name>(check_adadelta_convergence,
                                                      batch_ligands, population_number,
                                                      converged_ligands_b, ls_selected_b, scores_b,
                                                      best_scores_b, best_scores_ages_b,
                                                      chromosomes_b, best_chromosomes_b,
                                                      completed_iterations, finalize);
  }
} // namespace mudock
