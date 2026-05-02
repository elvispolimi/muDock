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
                                    const fp_type epsilon,
                                    const fp_type rho,
                                    gradient *__restrict__ gradients_b,
                                    chromosome *__restrict__ population_b,
                                    chromosome *__restrict__ adadelta_e_g2_b,
                                    chromosome *__restrict__ adadelta_e_dw2_b) {
    
    // Gradient size matches chromosome size (6 + max_rotamers)
    // TODO L is this ok?
    constexpr int gradient_size = 6 + max_static_bonds();
    const size_t total_individuals = batch_ligands * individuals_per_ligand;
    
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      gradient *gradients_l = gradients_b + ligand_index * individuals_per_ligand;
      chromosome *population_l = population_b + ligand_index * individuals_per_ligand;
      chromosome *e_g2_l = adadelta_e_g2_b + ligand_index * individuals_per_ligand;
      chromosome *e_dw2_l = adadelta_e_dw2_b + ligand_index * individuals_per_ligand;
      
      for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {
        gradient &grad = gradients_l[individual_index];
        chromosome &w = population_l[individual_index];
        chromosome &E_g2_i = e_g2_l[individual_index];
        chromosome &E_dw2_i = e_dw2_l[individual_index];
        
        // Apply AdaDelta update for each dimension
        for (int d = 0; d < gradient_size; ++d) {
          // Update E[g^2] running average: E[g^2] = rho * E[g^2] + (1-rho) * g^2
          E_g2_i[d] = rho * E_g2_i[d] + (1.0f - rho) * grad[d] * grad[d];
          
          // Compute RMS of gradient: RMS[g] = sqrt(E[g^2] + epsilon)
          const fp_type rms_g = std::sqrt(E_g2_i[d] + epsilon);
          
          // Compute RMS of delta (from previous step): RMS[delta] = sqrt(E[delta^2] + epsilon)
          const fp_type rms_dw = std::sqrt(E_dw2_i[d] + epsilon);
          
          // Compute delta_w: delta_w = -RMS[delta] / RMS[g] * g
          const fp_type delta_w = -(rms_dw / rms_g) * grad[d];
          
          // Update E[delta^2] running average: E[delta^2] = rho * E[delta^2] + (1-rho) * delta_w^2
          E_dw2_i[d] = rho * E_dw2_i[d] + (1.0f - rho) * delta_w * delta_w;
          
          // Update weights: w = w + delta_w
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
                                                ADADELTA_EPSILON,
                                                ADADELTA_RHO,
                                                gradients_b,
                                                population_b,
                                                adadelta_e_g2_b,
                                                adadelta_e_dw2_b);
  }

  // TODO L what to do with this? i moved the iterations in adadelta.hpp
  // template<>
  // void adadelta_kernel<queue_cpp>::operator()() {
  //
  // }
} // namespace mudock