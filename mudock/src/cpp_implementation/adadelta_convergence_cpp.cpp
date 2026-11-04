#include <mudock/type_alias.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>
#include <mudock/compute/adadelta_convergence.hpp>
#include <mudock/cpp_implementation/adadelta_convergence_cpp.hpp>
#include <cmath>

namespace mudock {
  void check_adadelta_convergence(const int batch_ligands,
                                  // const bool autostop,
                                  const int individuals_per_ligand,
                                  const int* __restrict__ converged_ligands_b,
                                  int* __restrict__ ls_selected_b,
                                  const fp_type *__restrict__ scores_b,
                                  fp_type* __restrict__ best_scores_b, 
                                  int* __restrict__ best_scores_ages_b) {
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      if (converged_ligands_b[ligand_index]) continue;

      const fp_type* __restrict__ scores_l     = scores_b + ligand_index * individuals_per_ligand;
      fp_type *__restrict__ best_scores_l      = best_scores_b + ligand_index * individuals_per_ligand;
      int* __restrict__ best_scores_ages_l = best_scores_ages_b + ligand_index * individuals_per_ligand;
      int* __restrict__ ls_selected_l      = ls_selected_b + ligand_index * individuals_per_ligand;

      for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {  
        if (!ls_selected_l[individual_index]) continue;
        
        if (scores_l[individual_index] < best_scores_l[individual_index]) {
          best_scores_l[individual_index] = scores_l[individual_index];
          best_scores_ages_l[individual_index] = 0;
        } else {
          best_scores_ages_l[individual_index]++;
        }

        if (best_scores_ages_l[individual_index] > MAX_BEST_SCORE_AGE) { // score didn't improve for n iterations 
          ls_selected_l[individual_index] = 0;
        }
      }

    }

  }

  
  template<>
  void adadelta_convergence_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<this->adadelta_convergence_region_name>(check_adadelta_convergence,
                                                            batch_ligands,
                                                            // autostop,
                                                            population_number,
                                                            converged_ligands_b,
                                                            ls_selected_b,
                                                            scores_b,
                                                            best_scores_b, 
                                                            best_scores_ages_b);
  }

} // namespace mudock