#include <mudock/type_alias.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/molecule/constraints.hpp>
#include <cmath>

// Gradient size matches chromosome size (6 + max_rotamers)
constexpr int gradient_size = 6 + max_static_bonds();

namespace mudock {
  // Inline AdaDelta update - applies the AdaDelta update rule to all individuals
  inline void apply_adadelta_update(const int batch_ligands,
                                    const int individuals_per_ligand,
                                    const fp_type epsilon,
                                    const fp_type rho,
                                    gradient *__restrict__ gradients_b,
                                    chromosome *__restrict__ population_b) {
    
    // AdaDelta state: running averages of gradient squared and delta squared
    // These need to persist across iterations - stored in thread-local storage
    // Format: E[g^2] and E[delta^2] for each dimension
    static thread_local std::vector<chromosome> E_g2;
    static thread_local std::vector<chromosome> E_dw2;
    static thread_local std::vector<int> initialized;
    
    const size_t total_individuals = batch_ligands * individuals_per_ligand;
    
    // Resize state buffers if needed
    if (E_g2.size() < total_individuals) {
      E_g2.resize(total_individuals);
      E_dw2.resize(total_individuals);
      initialized.resize(total_individuals, 0);
    }
    
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      gradient *gradients_l = gradients_b + ligand_index * individuals_per_ligand;
      chromosome *population_l = population_b + ligand_index * individuals_per_ligand;
      
      for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {
        const int idx = ligand_index * individuals_per_ligand + individual_index;
        
        gradient &grad = gradients_l[individual_index];
        chromosome &w = population_l[individual_index];
        
        chromosome &E_g2_i = E_g2[idx];
        chromosome &E_dw2_i = E_dw2[idx];
        
        // Initialize state if needed
        if (!initialized[idx]) {
          for (int d = 0; d < gradient_size; ++d) {
            E_g2_i[d] = 0.0f;
            E_dw2_i[d] = 0.0f;
          }
          initialized[idx] = 1;
        }
        
        // Apply AdaDelta update for each dimension
        for (int dim = 0; dim < gradient_size; ++dim) {
          // Update E[g^2] running average: E[g^2] = rho * E[g^2] + (1-rho) * g^2
          E_g2_i[dim] = rho * E_g2_i[dim] + (1.0f - rho) * grad[dim] * grad[dim];
          
          // Compute RMS of gradient: RMS[g] = sqrt(E[g^2] + epsilon)
          const fp_type rms_g = std::sqrt(E_g2_i[dim] + epsilon);
          
          // Compute RMS of delta (from previous step): RMS[delta] = sqrt(E[delta^2] + epsilon)
          const fp_type rms_dw = std::sqrt(E_dw2_i[dim] + epsilon);
          
          // Compute delta_w: delta_w = -RMS[delta] / RMS[g] * g
          const fp_type delta_w = -(rms_dw / rms_g) * grad[dim];
          
          // Update E[delta^2] running average: E[delta^2] = rho * E[delta^2] + (1-rho) * delta_w^2
          E_dw2_i[dim] = rho * E_dw2_i[dim] + (1.0f - rho) * delta_w * delta_w;
          
          // Update weights: w = w + delta_w
          w[dim] = w[dim] + delta_w;
        }
      }
    }
  }

  template<>
  void adadelta_kernel<queue_cpp>::compute_gradients(){
    this->score_stage->compute_gradient();
  }
  
  template<>
  void adadelta_kernel<queue_cpp>::apply_adadelta() {
    q->invoke_kernel<this->adadelta_region_name>(apply_adadelta_update,
                                                batch_ligands,
                                                individuals_per_ligand,
                                                ADADELTA_EPSILON,
                                                ADADELTA_RHO,
                                                gradients_b,
                                                population_b);
  }

  // TODO L what to do with this? i moved the iterations in adadelta.hpp
  // template<>
  // void adadelta_kernel<queue_cpp>::operator()() {
  //
  // }
} // namespace mudock