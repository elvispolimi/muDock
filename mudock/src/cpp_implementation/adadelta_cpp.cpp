#include <mudock/type_alias.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>

#define MAX_AD_ITERATIONS 300

namespace mudock {

  // TODO L add an early stop criterion based on convergence
  inline void perform_local_search_adadelta(const fp_type epsilon,
                                            const fp_type rho) {
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {
        
        // Adadelta algorithm
        fp_type num = 0;
        fp_type den = 0;
        chromosome w            = {0};
        chromosome delta_w      = {0};
        chromosome grad         = {0};
        chromosome E_gt         = {0};
        chromosome E_gt_past    = {0};
        chromosome E_delta_wt_past         = {0};
        chromosome E_delta_wt_past_past    = {0};
        int        num_dim      = sizeof(grad) / sizeof(grad[0]);

        // Loop over time steps
        for (int t = 0; t < MAX_AD_ITERATIONS; ++t) {
        
          // Loop over each dimension
          for (int dim = 0; dim < num_dim; ++dim) {

            E_delta_wt_past[dim] = rho*E_delta_wt_past_past[dim] + (1-rho)*square(delta_w[dim]);
            num = sqrt(E_delta_wt_past[dim]+epsilon);
  
            E_gt[dim] = rho*E_gt_past[dim] + (1-rho)*square(grad[dim]);
            den = sqrt(E_gt[dim]+epsilon);
            
            // Find delta_w
            delta_w[dim] = -(num / den) * grad[dim];
            
            // Update weights
            w[dim] = w[dim] + delta_w[dim];
          }
        }
      }
    }
  };

  template<>
  void adadelta_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<this->adt_region_name>(perform_local_search_adadelta
                                            /* other params */);
  }
} // namespace mudock
