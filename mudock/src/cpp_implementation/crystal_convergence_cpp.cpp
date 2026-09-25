#include <mudock/type_alias.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>
#include <mudock/compute/crystal_convergence.hpp>
#include <mudock/cpp_implementation/crystal_convergence_cpp.hpp>
#include <cmath>

namespace mudock {
  void check_convergence(const int batch_ligands,
                         const int batch_atoms,
                         const int individuals_per_ligand,
                         const int* __restrict__ num_atoms_b,
                         int* generation,
                         int* __restrict__ converged_ligands_b,
                         const fp_type* __restrict__ template_x_b,
                         const fp_type* __restrict__ template_y_b,
                         const fp_type* __restrict__ template_z_b,
                         const fp_type *__restrict__ x_scratch_b,
                         const fp_type *__restrict__ y_scratch_b,
                         const fp_type *__restrict__ z_scratch_b,
                         fp_type *__restrict__ rmsd_best_pose_b,
                         const fp_type *__restrict__ scores_b) {
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      const int converged_ligand = converged_ligands_b[ligand_index];
      if (converged_ligand) {
        continue;
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////
      // WARNING: CRITERION 1 AND 2 ARE ARTIFICAL CONVERGENCE CRITERION, NOT REAL ONE, BECAUSE THEY ASSUME
      //          WE HAVE THE CRYSTAL, WHICH IS OK FOR EVALUATING PERFORMANCE, BUT NOT USEFUL FOR REAL CASE
      //          SCENARIOS WHERE INSTEAD WE MUST CHECK CONVERGENCE LOOKING AT THE SCORE OR POPULATION
      //          STAGNATION ETC...
      //////////////////////////////////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////////////////////////////////
      // CRITERION 1: best_score <= (crystal_score + 1)
      //////////////////////////////////////////////////////////////////////////////////////////////////////
      // TODO L: implement 
      const fp_type *__restrict__ scores_l = scores_b + ligand_index * individuals_per_ligand;
      int best_index = 0;
      fp_type best_score = scores_l[best_index];
      for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {
        if (scores_l[individual_index] < best_score){
          best_index = individual_index;
          best_score = scores_l[individual_index];
        }
      }
      printf("Gen %d --- Best: %f\n", *generation, double(best_score));
      bool good_score = best_score <= fp_type{-15.0} ? true : false;


      //////////////////////////////////////////////////////////////////////////////////////////////////////
      // CRITERION 2: RMSD best scoring pose < 1 Ångström
      //////////////////////////////////////////////////////////////////////////////////////////////////////
      const int num_atoms   = num_atoms_b[ligand_index];
      const int atom_stride = ligand_index * batch_atoms;
      const fp_type* template_x_l = template_x_b + atom_stride;
      const fp_type* template_y_l = template_y_b + atom_stride;
      const fp_type* template_z_l = template_z_b + atom_stride;

      const fp_type *__restrict__ scratch_x = x_scratch_b + atom_stride * individuals_per_ligand;
      const fp_type *__restrict__ scratch_y = y_scratch_b + atom_stride * individuals_per_ligand;
      const fp_type *__restrict__ scratch_z = z_scratch_b + atom_stride * individuals_per_ligand;
      const fp_type *__restrict__ scratch_x_l = scratch_x + best_index * batch_atoms;
      const fp_type *__restrict__ scratch_y_l = scratch_y + best_index * batch_atoms;
      const fp_type *__restrict__ scratch_z_l = scratch_z + best_index * batch_atoms;

      rmsd_best_pose_b[ligand_index] = compute_rmsd(template_x_l, template_y_l, template_z_l,
                                                      scratch_x_l, scratch_y_l, scratch_z_l, 
                                                      num_atoms);
      // fp_type rmsd_best_scoring_pose = 0;
      bool good_rmsd = rmsd_best_pose_b[ligand_index] < 1 ? true : false;

      //////////////////////////////////////////////////////////////////////////////////////////////////////
      // mark convergence
      //////////////////////////////////////////////////////////////////////////////////////////////////////
      if (good_score || good_rmsd) {
        converged_ligands_b[ligand_index] = *generation;
        printf("CRYSTAL FOUND! score = %f, rmsd = %f\n", best_score, rmsd_best_pose_b[ligand_index]);
      }

    }
    *generation += 1;
  }

  
  template<>
  void crystal_convergence_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<this->crystal_convergence_region_name>(check_convergence,
                                                            batch_ligands,
                                                            batch_atoms,
                                                            population_number,
                                                            num_atoms_b,
                                                            &current_generation,
                                                            converged_ligands_b,
                                                            template_x_b,
                                                            template_y_b,
                                                            template_z_b,
                                                            x_scratch_b,
                                                            y_scratch_b,
                                                            z_scratch_b,
                                                            rmsd_best_pose_b,
                                                            scores_b);
  }

} // namespace mudock