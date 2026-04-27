#include <mudock/type_alias.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/compute/adt_score_kernel.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/molecule/constraints.hpp>
#include <cmath>

// AdaDelta hyperparameters
#define ADADELTA_RHO 0.95f
#define ADADELTA_EPSILON 1e-6f
#define MAX_AD_ITERATIONS 300

// Gradient size matches chromosome size (6 + max_rotamers)
constexpr int gradient_size = 6 + max_static_bonds();

namespace mudock {

  /*
  // Inline gradient computation - reuses the same logic as adt_score
  // Computes gradients for ALL individuals in the batch in one call
  // This is a simplified version - for full accuracy, reuse calc_gradient from adt_score_cpp.cpp
  inline void compute_adadelta_gradients(const int batch_atoms,
                                         const int batch_ligands,
                                         const int individuals_per_ligand,
                                         const fp_type *__restrict__ x_scratch_b,
                                         const fp_type *__restrict__ y_scratch_b,
                                         const fp_type *__restrict__ z_scratch_b,
                                         const fp_type *__restrict__ vols_b,
                                         const fp_type *__restrict__ solpars_b,
                                         const fp_type *__restrict__ charges_b,
                                         const int *__restrict__ num_atoms_b,
                                         const int *__restrict__ num_rotamers_b,
                                         const int *__restrict__ num_nonbonds_b,
                                         const int *__restrict__ nonbond_a1_b,
                                         const int *__restrict__ nonbond_a2_b,
                                         const fp_type *__restrict__ nonbond_cA_b,
                                         const fp_type *__restrict__ nonbond_cB_b,
                                         const int *__restrict__ nonbond_xB_b,
                                         const fp_type *__restrict__ grid_maps,
                                         const fp_type *__restrict__ minimum,
                                         const fp_type *__restrict__ maximum,
                                         const fp_type *__restrict__ center,
                                         const int *__restrict__ map_offsets_b,
                                         const int map_index_x,
                                         const int map_index_xy,
                                         const int map_index_xyz,
                                         gradient *__restrict__ gradients_b) {
    // This function computes gradients for ALL individuals in the batch
    // The gradient structure matches the chromosome structure (6 + num_rotamers dimensions)
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      const int atom_stride  = ligand_index * batch_atoms;
      const int num_atoms    = num_atoms_b[ligand_index];
      const int num_nonbonds = num_nonbonds_b[ligand_index + 1] - num_nonbonds_b[ligand_index];
      const int num_rotamers = num_rotamers_b[ligand_index];

      const fp_type *__restrict__ scratch_x = x_scratch_b + atom_stride * individuals_per_ligand;
      const fp_type *__restrict__ scratch_y = y_scratch_b + atom_stride * individuals_per_ligand;
      const fp_type *__restrict__ scratch_z = z_scratch_b + atom_stride * individuals_per_ligand;
      const fp_type *__restrict__ vol_l     = vols_b + atom_stride;
      const fp_type *__restrict__ solpar_l  = solpars_b + atom_stride;
      const fp_type *__restrict__ charge_l  = charges_b + atom_stride;
      const int *__restrict__ map_offsets_l = map_offsets_b + atom_stride;
      const int *__restrict__ nonbond_a1_l  = nonbond_a1_b + num_nonbonds_b[ligand_index];
      const int *__restrict__ nonbond_a2_l  = nonbond_a2_b + num_nonbonds_b[ligand_index];
      const fp_type *nonbond_cA_l           = nonbond_cA_b + num_nonbonds_b[ligand_index];
      const fp_type *nonbond_cB_l           = nonbond_cB_b + num_nonbonds_b[ligand_index];
      const int *nonbond_xB_l               = nonbond_xB_b + num_nonbonds_b[ligand_index];

      gradient *__restrict__ gradients_l = gradients_b + ligand_index * individuals_per_ligand;
      
      // Compute gradient for each individual in the population
      for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {
        const fp_type *__restrict__ scratch_x_l = scratch_x + individual_index * batch_atoms;
        const fp_type *__restrict__ scratch_y_l = scratch_y + individual_index * batch_atoms;
        const fp_type *__restrict__ scratch_z_l = scratch_z + individual_index * batch_atoms;
        
        gradient &grad = gradients_l[individual_index];
        
        // Initialize gradient to zero
        for (int d = 0; d < gradient_size; ++d) {
          grad[d] = 0.0f;
        }

        // Compute gradient for each atom using trilinear interpolation
        // This replicates the logic from calc_gradient in adt_score_cpp.cpp
        for (int atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
          fp_type coord[3]{scratch_x_l[atom_idx], scratch_y_l[atom_idx], scratch_z_l[atom_idx]};

          const auto diff_x = coord[0] - center[0];
          const auto diff_y = coord[1] - center[1];
          const auto diff_z = coord[2] - center[2];

          // Check if atom is within grid bounds
          if (coord[0] < minimum[0] || coord[0] > maximum[0] || coord[1] < minimum[1] ||
              coord[1] > maximum[1] || coord[2] < minimum[2] || coord[2] > maximum[2]) {
            // Outside grid - simple penalty gradient
            const fp_type dist = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
            const fp_type epenalty = dist * ENERGYPENALTY;
            // Gradient points toward center
            grad[0] += 2.0f * diff_x * epenalty;  // dE/dx
            grad[1] += 2.0f * diff_y * epenalty;  // dE/dy
            grad[2] += 2.0f * diff_z * epenalty;  // dE/dz
          } else {
            // Inside grid - compute trilinear interpolation gradients
            const auto &atom_charge = charge_l[atom_idx];
            
            // Normalize coordinates to grid space
            fp_type gcoord[3];
            gcoord[0] = (coord[0] - minimum[0]) * inv_spacing;
            gcoord[1] = (coord[1] - minimum[1]) * inv_spacing;
            gcoord[2] = (coord[2] - minimum[2]) * inv_spacing;

            const int u0 = static_cast<int>(gcoord[0]);
            const int v0 = static_cast<int>(gcoord[1]);
            const int w0 = static_cast<int>(gcoord[2]);
            
            const fp_type pu[2] = {fp_type{1} - (gcoord[0] - u0), gcoord[0] - u0};
            const fp_type pv[2] = {fp_type{1} - (gcoord[1] - v0), gcoord[1] - v0};
            const fp_type pw[2] = {fp_type{1} - (gcoord[2] - w0), gcoord[2] - w0};

            // Get map indices
            const int map_offset = map_offsets_l[atom_idx];
            const fp_type *elec_map = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
            const fp_type *desolv_map = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);
            const fp_type *atom_map = grid_maps + map_offset;

            // Compute gradient contributions using trilinear interpolation
            // dE/dx = sum over 8 corners of (dE/dgrid * dgrid/dx)
            // dgrid/dx = weight derivative = +/- inv_spacing
            for (int i = 0; i < 2; ++i) {
              for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                  const int idx = (u0 + i) + (v0 + j) * map_index_x + (w0 + k) * map_index_xy;
                  const fp_type weight = pu[i] * pv[j] * pw[k];
                  const fp_type weight_deriv_x = (i == 0 ? -1.0f : 1.0f) * inv_spacing * pv[j] * pw[k];
                  const fp_type weight_deriv_y = (j == 0 ? -1.0f : 1.0f) * inv_spacing * pu[i] * pw[k];
                  const fp_type weight_deriv_z = (k == 0 ? -1.0f : 1.0f) * inv_spacing * pu[i] * pv[j];
                  
                  // Electrostatic gradient contribution
                  const fp_type eval = (idx < map_index_xyz) ? elec_map[idx] : 0.0f;
                  grad[0] += weight_deriv_x * eval * atom_charge;
                  grad[1] += weight_deriv_y * eval * atom_charge;
                  grad[2] += weight_deriv_z * eval * atom_charge;
                  
                  // Energy map gradient contribution
                  const fp_type eval_emap = (idx < map_index_xyz) ? atom_map[idx] : 0.0f;
                  grad[0] += weight_deriv_x * eval_emap;
                  grad[1] += weight_deriv_y * eval_emap;
                  grad[2] += weight_deriv_z * eval_emap;
                  
                  // Desolvation gradient contribution
                  const fp_type dval = (idx < map_index_xyz) ? desolv_map[idx] : 0.0f;
                  const fp_type vol = vol_l[atom_idx];
                  const fp_type solp = solpar_l[atom_idx];
                  grad[0] += weight_deriv_x * vol * (dval + solp * std::fabs(atom_charge));
                  grad[1] += weight_deriv_y * vol * (dval + solp * std::fabs(atom_charge));
                  grad[2] += weight_deriv_z * vol * (dval + solp * std::fabs(atom_charge));
                }
              }
            }
          }
        }

        // Add internal (non-bonded) energy gradients
        // This is a simplified version - full implementation would replicate the logic from calc_gradient
        if (num_rotamers > 0) {
          for (int i = 0; i < num_nonbonds; ++i) {
            const int &a1 = nonbond_a1_l[i];
            const int &a2 = nonbond_a2_l[i];

            const auto diff_x = scratch_x_l[a1] - scratch_x_l[a2];
            const auto diff_y = scratch_y_l[a1] - scratch_y_l[a2];
            const auto diff_z = scratch_z_l[a1] - scratch_z_l[a2];
            const fp_type distance_two = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
            const fp_type distance = std::sqrt(std::max(distance_two, RMIN_ELEC_SQUARE));
            
            if (distance_two < nbc2) {
              // Simplified VdW gradient (full version would use lookup tables)
              const fp_type inv_r = 1.0f / distance;
              const fp_type dir_x = diff_x * inv_r;
              const fp_type dir_y = diff_y * inv_r;
              const fp_type dir_z = diff_z * inv_r;
              
              // Placeholder for VdW gradient - would need proper coefficient lookup
              const fp_type dE_dr_vdw = 0.01f * distance;  // Simplified
              
              grad[0] += dE_dr_vdw * dir_x;
              grad[1] += dE_dr_vdw * dir_y;
              grad[2] += dE_dr_vdw * dir_z;
            }
          }
        }

        // Map atom position gradients to chromosome gene gradients
        // Translation genes (0-2): already computed as sum of atom gradients
        // Rotation genes (3-5): require Jacobian computation (simplified here)
        // For now, translation gradients are complete, rotation placeholders
        for (int rot_idx = 0; rot_idx < num_rotamers; ++rot_idx) {
          grad[6 + rot_idx] = 0.001f;  // Placeholder - requires full rotational Jacobian
        }
      }
    }
  }*/

  // template<>
  // void adadelta_kernel<queue_cpp>::compute_gradients() {
  //   q->invoke_kernel<this->gradient_region_name>(/*compute_adadelta_gradients,*/
  //                                                batch_atoms,
  //                                                batch_ligands,
  //                                                scores_per_ligand,
  //                                                x_scratch_b,
  //                                                y_scratch_b,
  //                                                z_scratch_b,
  //                                                vols_b,
  //                                                solpars_b,
  //                                                charges_b,
  //                                                num_atoms_b,
  //                                                num_rotamers_b,
  //                                                num_nonbonds_b,
  //                                                nonbond_a1_b,
  //                                                nonbond_a2_b,
  //                                                nonbond_cA_b,
  //                                                nonbond_cB_b,
  //                                                nonbond_xB_b,
  //                                                grid_maps,
  //                                                minimum,
  //                                                maximum,
  //                                                center,
  //                                                map_offsets_b,
  //                                                map_index_x,
  //                                                map_index_xy,
  //                                                map_index_xyz,
  //                                                gradients_b);
  // }

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
    // Apply the AdaDelta update using the gradients and population
    // population_b should be set via set_population() before calling this
    q->invoke_kernel<this->adt_region_name>(apply_adadelta_update,
                                            batch_ligands,
                                            scores_per_ligand,
                                            ADADELTA_EPSILON,
                                            ADADELTA_RHO,
                                            gradients_b,
                                            population_b);
  }

  // TODO L what to do with this? i moved the iterations in adadelta.hpp
  // template<>
  // void adadelta_kernel<queue_cpp>::operator()() {
  //   for (int i = 0; i < MAX_AD_ITERATIONS; ++i){
  //     compute_gradients();
  //     apply_adadelta();
  //   }
  // }
} // namespace mudock