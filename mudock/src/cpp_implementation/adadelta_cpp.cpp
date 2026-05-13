#include <mudock/type_alias.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>
#include <mudock/compute/adadelta_kernel.hpp>
#include <mudock/compute/adadelta.hpp>
#include <cmath>

namespace mudock {
  // Inline AdaDelta update - applies the AdaDelta update rule to all individuals
  inline void apply_adadelta_update(const int batch_ligands,
                                    const int individuals_per_ligand,
                                    gradient *__restrict__ gradients_b,
                                    chromosome *__restrict__ population_b,
                                    chromosome *__restrict__ adadelta_e_g2_b,
                                    chromosome *__restrict__ adadelta_e_dw2_b,
                                    int *__restrict__ stall_counter_b,
                                    int *__restrict__ inactive_b,
                                    int i) { // TODO L remove i (number iteration) if not needed 
    
    // Gradient size matches chromosome size (6 + max_rotamers)
    // TODO L is this ok?
    constexpr int gradient_size = 6 + max_static_bonds();
    
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      gradient *gradients_l = gradients_b + ligand_index * individuals_per_ligand;
      chromosome *population_l = population_b + ligand_index * individuals_per_ligand;
      chromosome *e_g2_l = adadelta_e_g2_b + ligand_index * individuals_per_ligand;
      chromosome *e_dw2_l = adadelta_e_dw2_b + ligand_index * individuals_per_ligand;
      int *stall_counter_l = stall_counter_b + ligand_index * individuals_per_ligand;
      int *inactive_l = inactive_b + ligand_index * individuals_per_ligand;
      
      for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {
        gradient &grad = gradients_l[individual_index];
        chromosome &w = population_l[individual_index];
        chromosome &E_g2_i = e_g2_l[individual_index];
        chromosome &E_dw2_i = e_dw2_l[individual_index];
        int &stall_counter = stall_counter_l[individual_index];
        int &inactive = inactive_l[individual_index];

        // Reset values for in each generation of genetic.
        // TODO L this reset is crap. make it better
        if (i == 0){
          stall_counter = 0;
          inactive = 0;
        }
        
        // Skip individual if it already converged 
        if (inactive == 1){
          //printf("id: %d iter: %d, skipping for ES...\n", individual_index, i);
          continue;
        }

        fp_type delta_norm2 = 0;

        // Apply AdaDelta update for each dimension
        for (int d = 0; d < gradient_size; ++d) {
          E_g2_i[d] = ADADELTA_RHO * E_g2_i[d] + (1.0f - ADADELTA_RHO) * grad[d] * grad[d];
          
          const fp_type rms_g = std::sqrt(E_g2_i[d] + ADADELTA_EPSILON);
          
          const fp_type rms_dw = std::sqrt(E_dw2_i[d] + ADADELTA_EPSILON);
          
          const fp_type delta_w = -(rms_dw / rms_g) * grad[d];

          delta_norm2 += delta_w * delta_w;
          
          E_dw2_i[d] = ADADELTA_RHO * E_dw2_i[d] + (1.0f - ADADELTA_RHO) * delta_w * delta_w;

          w[d] = w[d] + delta_w;
        }

        // per la convergenza: calcolo l'uno per cento dello spazio di movimento delle varie dimensioni (-9;9 A per traslazioni, -180;180 per gli angoli).
        // controllo se esiste almeno una delle dimensioni che ha un cambiamento maggiore di questo. converge quando tutte le dimensioni sono sotto la loro soglia

        // // todo 
        // - usare address sanitizer
        // - gradiente della dimensione giusta e non massima
        // - controllare se è VdW a causare l'esplosione dello score dopo tot iterazioni
        // - chiedere a gianmarco se ha gia fatto to_mol2

        if(i % 5 == 0){
          printf("id: %d iter: %d, delta_norm2: %f\n", individual_index, i, delta_norm2);
        }

        // Checking convergence
        if (delta_norm2 < ADADELTA_CONVERGENCE_THRESHOLD) {
          stall_counter++;
        }
        else{
          stall_counter = 0;
        }
        if (stall_counter >= ADADELTA_CONVERGENCE_PATIENCE){
          inactive = 1;
          printf("HIT CONVERGENCE - id: %d at iter: %d\n", individual_index, i);
        }

      }
    }
  }

  template<>
  void adadelta_kernel<queue_cpp>::compute_gradients(){
    (this->score_stage).get()->compute_gradient();
  }
  
  template<>
  void adadelta_kernel<queue_cpp>::apply_adadelta(int i) {
    q->invoke_kernel<this->adadelta_region_name>(apply_adadelta_update,
                                                batch_ligands,
                                                individuals_per_ligand,
                                                gradients_b,
                                                population_b,
                                                adadelta_e_g2_b,
                                                adadelta_e_dw2_b,
                                                stall_counter_b,
                                                inactive_b,
                                                i);
  }

  // TODO L what to do with this? i moved the iterations in adadelta.hpp
  // template<>
  // void adadelta_kernel<queue_cpp>::operator()() {
  //
  // }
} // namespace mudock