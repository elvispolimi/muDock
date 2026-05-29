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
                                    int *__restrict__ active_individuals_b) {
    
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      chromosome *__restrict__ population_l         = population_b + ligand_index * individuals_per_ligand;
      chromosome *__restrict__ e_g2_l               = adadelta_e_g2_b + ligand_index * individuals_per_ligand;
      chromosome *__restrict__ e_dw2_l              = adadelta_e_dw2_b + ligand_index * individuals_per_ligand;
      gradient   *__restrict__ gradients_l          = gradients_b + ligand_index * individuals_per_ligand;
      const int  *__restrict__ active_individuals_l = active_individuals_b + ligand_index * individuals_per_ligand;
      
      const int num_rotamers = num_rotamers_b[ligand_index];
      const int gradient_size = 6 + num_rotamers;
      
      for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {
        if (!active_individuals_l[individual_index]){
          continue;
        }
        
        gradient   &grad = gradients_l[individual_index];
        chromosome &w = population_l[individual_index];
        chromosome &E_g2_i = e_g2_l[individual_index];
        chromosome &E_dw2_i = e_dw2_l[individual_index];

        // Apply AdaDelta update for each dimension
        for (int d = 0; d < gradient_size; ++d) {
          E_g2_i[d] = RHO * E_g2_i[d] + (1.0f - RHO) * grad[d] * grad[d];
          
          const fp_type rms_g = std::sqrt(E_g2_i[d] + EPSILON);
          
          const fp_type rms_dw = std::sqrt(E_dw2_i[d] + EPSILON);
          
          fp_type delta_w = -(rms_dw / rms_g) * grad[d];

          // TODO L understand if we should use clamping
          
          // delta_w = std::clamp(delta_w, -MAX_STEP, MAX_STEP);
          
          // if (d < 3)      delta_w = std::clamp(delta_w, -MAX_STEP_POS, MAX_STEP_POS);
          // else if (d < 6) delta_w = std::clamp(delta_w, -MAX_STEP_ROT, MAX_STEP_ROT);
          // else            delta_w = std::clamp(delta_w, -MAX_STEP_TORS, MAX_STEP_TORS);
          
          E_dw2_i[d] = RHO * E_dw2_i[d] + (1.0f - RHO) * delta_w * delta_w;

          w[d] = w[d] + delta_w;

        }
      }
    }
  }

  template<>
  void adadelta_kernel<queue_cpp>::compute_gradients(){
    (this->score_stage).get()->compute_gradient();
  }
  
  template<>
  void adadelta_kernel<queue_cpp>::apply_adadelta() {
    q->invoke_kernel<this->adadelta_region_name>(apply_adadelta_update,
                                                batch_ligands,
                                                individuals_per_ligand,
                                                gradients_b,
                                                population_b,
                                                num_rotamers_b,
                                                adadelta_e_g2_b,
                                                adadelta_e_dw2_b,
                                                active_individuals_b);
  }

  // TODO L what to do with this? i moved the iterations in adadelta.hpp
  // template<>
  // void adadelta_kernel<queue_cpp>::operator()() {
  //
  // }
} // namespace mudock