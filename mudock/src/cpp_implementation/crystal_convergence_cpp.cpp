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
                                fp_type *__restrict__ scores_b) {
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      
      fp_type *__restrict__ scores_l = scores_b + ligand_index * individuals_per_ligand;
      
      fp_type best_score = scores_l[0];
      for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {
        if (scores_l[individual_index] < best_score){
          best_score = scores_l[individual_index];
        }
      }
      
      if (best_score <= fp_type{-15.0}) {
        printf("CRYSTAL FOUND!\n");
      }

    }
  }

  
  template<>
  void crystal_convergence_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<this->crystal_convergence_region_name>(check_convergence,
                                                            batch_ligands,
                                                            population_number,
                                                            scores_b);
  }

} // namespace mudock