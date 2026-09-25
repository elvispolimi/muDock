#include <mudock/type_alias.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>
#include <mudock/compute/crystal_convergence.hpp>
#include <mudock/cpp_implementation/crystal_convergence_cpp.hpp>
#include <cmath>

namespace mudock {
  inline void check_convergence(const int batch_ligands,
                                const int individuals_per_ligand,
                                int* generation,
                                int* __restrict__ converged_ligands_b,
                                fp_type *__restrict__ scores_b) {
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      const int converged_ligand = converged_ligands_b[ligand_index];
      if (converged_ligand) {
        continue;
      }
      fp_type *__restrict__ scores_l = scores_b + ligand_index * individuals_per_ligand;
      
      fp_type best_score = scores_l[0];
      for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {
        if (scores_l[individual_index] < best_score){
          best_score = scores_l[individual_index];
        }
      }
      
      // if (autostop && has_converged(ligand_index, for_how_long_best_b, tolerance_window, best_so_far_b[ligand_index], best, best_score_diff_thld, scores, population_number, score_variance_thld)) {
      //   converged_ligands_b[ligand_index] = generation;
      // }

      if (best_score <= fp_type{-15.0}) {
        converged_ligands_b[ligand_index] = *generation;
        printf("CRYSTAL FOUND!\n");
      }

    }
    *generation += 1;
  }

  
  template<>
  void crystal_convergence_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<this->crystal_convergence_region_name>(check_convergence,
                                                            batch_ligands,
                                                            population_number,
                                                            &current_generation,
                                                            converged_ligands_b,
                                                            scores_b);
  }

} // namespace mudock