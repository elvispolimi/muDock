#include <mudock/type_alias.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>
#include <mudock/compute/adadelta_kernel.hpp>
#include <mudock/compute/adadelta.hpp>
#include <mudock/cpp_implementation/adadelta_cpp.hpp>
#include <cmath>

namespace mudock {
  // Inline AdaDelta update - applies the AdaDelta update rule to all individuals
  inline void apply_adadelta_update(const int batch_ligands,
                                    const int individuals_per_ligand,
                                    gradient *__restrict__ gradients_b,
                                    chromosome *__restrict__ population_b,
                                    int* __restrict__ num_rotamers_b,
                                    chromosome *__restrict__ adadelta_e_g2_b,
                                    chromosome *__restrict__ adadelta_e_dw2_b,
                                    int *__restrict__ stall_counter_b,
                                    int *__restrict__ inactive_b,
                                    int i,
                                    int convergence_patience,
                                    bool use_early_stopping) { // TODO L remove i (number iteration) if not needed 
    
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      chromosome *__restrict__ population_l    = population_b + ligand_index * individuals_per_ligand;
      chromosome *__restrict__ e_g2_l          = adadelta_e_g2_b + ligand_index * individuals_per_ligand;
      chromosome *__restrict__ e_dw2_l         = adadelta_e_dw2_b + ligand_index * individuals_per_ligand;
      gradient   *__restrict__ gradients_l     = gradients_b + ligand_index * individuals_per_ligand;
      int        *__restrict__ stall_counter_l = stall_counter_b + ligand_index * individuals_per_ligand;
      int        *__restrict__ inactive_l      = inactive_b + ligand_index * individuals_per_ligand;
      
      const int num_rotamers = num_rotamers_b[ligand_index];
      const int gradient_size = 6 + num_rotamers;
      
      for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {
        gradient   &grad = gradients_l[individual_index];
        chromosome &w = population_l[individual_index];
        chromosome &E_g2_i = e_g2_l[individual_index];
        chromosome &E_dw2_i = e_dw2_l[individual_index];
        int &stall_counter = stall_counter_l[individual_index];
        int &inactive = inactive_l[individual_index];

        // Reset values for in each generation of genetic.
        // TODO L this reset is garbage. make it better
        if (i == 0){
          stall_counter = 0;
          inactive = 0;
        }
        
        // Skip individual if it already converged 
        if (inactive == 1){
          //printf("id: %d iter: %d, skipping for ES...\n", individual_index, i);
          continue;
        }

        bool improvement = false;
      
        // Apply AdaDelta update for each dimension
        for (int d = 0; d < gradient_size; ++d) {
          E_g2_i[d] = RHO * E_g2_i[d] + (1.0f - RHO) * grad[d] * grad[d];
          
          const fp_type rms_g = std::sqrt(E_g2_i[d] + EPSILON);
          
          const fp_type rms_dw = std::sqrt(E_dw2_i[d] + EPSILON);
          
          fp_type delta_w = -(rms_dw / rms_g) * grad[d];

          // delta_w = std::clamp(delta_w, -MAX_STEP, MAX_STEP);
          
          if (d < 3)      delta_w = std::clamp(delta_w, -MAX_STEP_POS, MAX_STEP_POS);
          else if (d < 6) delta_w = std::clamp(delta_w, -MAX_STEP_ROT, MAX_STEP_ROT);
          else            delta_w = std::clamp(delta_w, -MAX_STEP_TORS, MAX_STEP_TORS);

          
          E_dw2_i[d] = RHO * E_dw2_i[d] + (1.0f - RHO) * delta_w * delta_w;

          w[d] = w[d] + delta_w;

          // To check convergence of an individual, we check if all the dimension has not improved significantly.
          // If at least one does improve, we continue the search.
          if (d < 3) {
            if (std::abs(delta_w) > CONVERGENCE_THRESHOLD_COORD) {
              improvement = true;
            }
          } else {
            if (std::abs(delta_w) > CONVERGENCE_THRESHOLD_ANGLE) {
              improvement = true;
            }
          }
        }

        // Checking convergence
        if (!improvement) {
          stall_counter++;
        }
        else{
          stall_counter = 0;
        }
        if (use_early_stopping && stall_counter >= convergence_patience){
          // TODO WARNING EARLY STOPPING IS DISABLED, UNCOMMENT INSTRUCTION TO ENABLE IT!
          inactive = 1;
          //printf("HIT CONVERGENCE - id: %d at iter: %d\n", individual_index, i);
          printf(".");
        }

      }
    }
  }

  template<>
  void adadelta_kernel<queue_cpp>::compute_gradients(){
    (this->score_stage).get()->compute_gradient();
  }
  
  template<>
  void adadelta_kernel<queue_cpp>::apply_adadelta(int i, int convergence_patience, bool use_early_stopping) {
    q->invoke_kernel<this->adadelta_region_name>(apply_adadelta_update,
                                                batch_ligands,
                                                individuals_per_ligand,
                                                gradients_b,
                                                population_b,
                                                num_rotamers_b,
                                                adadelta_e_g2_b,
                                                adadelta_e_dw2_b,
                                                stall_counter_b,
                                                inactive_b,
                                                i,
                                                convergence_patience,
                                                use_early_stopping);
  }

  // TODO L what to do with this? i moved the iterations in adadelta.hpp
  // template<>
  // void adadelta_kernel<queue_cpp>::operator()() {
  //
  // }
} // namespace mudock