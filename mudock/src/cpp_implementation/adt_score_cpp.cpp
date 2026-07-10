#include <cstring>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/grid/pi.hpp>
#include <mudock/cpp_implementation/adt_score_cpp.hpp>

#define FLATTENED_2D(x, y, index_x)              ((y) * index_x + (x))
#define FLATTENED_3D(x, y, z, index_x, index_xy) (index_xy * (z) + (y) * index_x + (x))

namespace mudock {
  namespace {
    fp_type trilinear_interpolation(const fp_type *__restrict__ map,
                                    const fp_type *__restrict__ coeffs,
                                    const int &map_index_x,
                                    const int &map_index_xy) {
    fp_type value{0};

    value = coeffs[0] * map[0] + value;
    value = coeffs[1] * map[map_index_xy] + value;
    value = coeffs[2] * map[map_index_x] + value;
    value = coeffs[3] * map[map_index_x + map_index_xy] + value;
    value = coeffs[4] * map[1] + value;
    value = coeffs[5] * map[1 + map_index_xy] + value;
    value = coeffs[6] * map[1 + map_index_x] + value;
    value = coeffs[7] * map[1 + map_index_x + map_index_xy] + value;

    return value;
    }

  inline void get_grid_values(const fp_type *__restrict__ map,
                              const int &map_index_x,
                              const int &map_index_xy,
                              fp_type *__restrict__ out_values) {
    out_values[0] = map[0];
    out_values[1] = map[map_index_xy];
    out_values[2] = map[map_index_x];
    out_values[3] = map[map_index_x + map_index_xy];
    out_values[4] = map[1];
    out_values[5] = map[1 + map_index_xy];
    out_values[6] = map[1 + map_index_x];
    out_values[7] = map[1 + map_index_x + map_index_xy];
  }

  void calc_energy(const int batch_atoms,
                   const int batch_ligands,
                   const int scores_per_ligand,
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
                   fp_type *__restrict__ scores_b) {
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      const int atom_stride  = ligand_index * batch_atoms;
      const int num_atoms    = num_atoms_b[ligand_index];
      const int num_nonbonds = num_nonbonds_b[ligand_index + 1] - num_nonbonds_b[ligand_index];
      const int num_rotamers = num_rotamers_b[ligand_index];

      const fp_type *__restrict__ scratch_x = x_scratch_b + atom_stride * scores_per_ligand;
      const fp_type *__restrict__ scratch_y = y_scratch_b + atom_stride * scores_per_ligand;
      const fp_type *__restrict__ scratch_z = z_scratch_b + atom_stride * scores_per_ligand;
      const fp_type *__restrict__ vol_l     = vols_b + atom_stride;
      const fp_type *__restrict__ solpar_l  = solpars_b + atom_stride;
      const fp_type *__restrict__ charge_l  = charges_b + atom_stride;
      const int *__restrict__ map_offsets_l = map_offsets_b + atom_stride;
      const int *__restrict__ nonbond_a1_l  = nonbond_a1_b + num_nonbonds_b[ligand_index];
      const int *__restrict__ nonbond_a2_l  = nonbond_a2_b + num_nonbonds_b[ligand_index];
      const fp_type *nonbond_cA_l           = nonbond_cA_b + num_nonbonds_b[ligand_index];
      const fp_type *nonbond_cB_l           = nonbond_cB_b + num_nonbonds_b[ligand_index];
      const int *nonbond_xB_l               = nonbond_xB_b + num_nonbonds_b[ligand_index];

      fp_type *__restrict__ scores_l = scores_b + ligand_index * scores_per_ligand;
      for (int scores_index = 0; scores_index < scores_per_ligand; ++scores_index) {
        const fp_type *__restrict__ scratch_x_l = scratch_x + scores_index * batch_atoms;
        const fp_type *__restrict__ scratch_y_l = scratch_y + scores_index * batch_atoms;
        const fp_type *__restrict__ scratch_z_l = scratch_z + scores_index * batch_atoms;

        fp_type elect_total_trilinear = 0;
        fp_type emap_total_trilinear  = 0;
        fp_type dmap_total_trilinear  = 0;
        const fp_type *electro_map = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
        const fp_type *desolv_map  = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

#pragma omp simd
        for (int index = 0; index < num_atoms; ++index) {
          fp_type coord[3]{scratch_x_l[index], scratch_y_l[index], scratch_z_l[index]};

          if (coord[0] < minimum[0] || coord[0] > maximum[0] || coord[1] < minimum[1] ||
              coord[1] > maximum[1] || coord[2] < minimum[2] || coord[2] > maximum[2]) {
            const auto diff_x      = coord[0] - center[0];
            const auto diff_y      = coord[1] - center[1];
            const auto diff_z      = coord[2] - center[2];
            const fp_type dist     = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
            const fp_type epenalty = dist * ENERGYPENALTY;
            elect_total_trilinear += epenalty;
            emap_total_trilinear += epenalty;
          } else {
            const auto &atom_charge = charge_l[index];
            const fp_type *atom_map = grid_maps + map_offsets_l[index];

            coord[0] = (coord[0] - minimum[0]) * inv_spacing;
            coord[1] = (coord[1] - minimum[1]) * inv_spacing;
            coord[2] = (coord[2] - minimum[2]) * inv_spacing;

            const int u0      = static_cast<int>(coord[0]);
            const fp_type p0u = coord[0] - static_cast<fp_type>(u0);
            const fp_type p1u = fp_type{1} - p0u;

            const int v0      = static_cast<int>(coord[1]);
            const fp_type p0v = coord[1] - static_cast<fp_type>(v0);
            const fp_type p1v = fp_type{1} - p0v;

            const int w0      = static_cast<int>(coord[2]);
            const fp_type p0w = coord[2] - static_cast<fp_type>(w0);
            const fp_type p1w = fp_type{1} - p0w;

            const fp_type pu[2] = {p1u, p0u};
            const fp_type pv[2] = {p1v, p0v};
            const fp_type pw[2] = {p1w, p0w};

            // Compute coefficients
            const fp_type coeffs[8] = {pu[0] * pv[0] * pw[0],
                                       pu[0] * pv[0] * pw[1],
                                       pu[0] * pv[1] * pw[0],
                                       pu[0] * pv[1] * pw[1],
                                       pu[1] * pv[0] * pw[0],
                                       pu[1] * pv[0] * pw[1],
                                       pu[1] * pv[1] * pw[0],
                                       pu[1] * pv[1] * pw[1]};

            // Precompute flattened indices
            const int base_index = FLATTENED_3D(u0, v0, w0, map_index_x, map_index_xy);
            // Trilinear Interpolationp
            elect_total_trilinear +=
                trilinear_interpolation(electro_map + base_index, coeffs, map_index_x, map_index_xy) *
                atom_charge;
            emap_total_trilinear +=
                trilinear_interpolation(atom_map + base_index, coeffs, map_index_x, map_index_xy);
            dmap_total_trilinear +=
                trilinear_interpolation(desolv_map + base_index, coeffs, map_index_x, map_index_xy) *
                std::fabs(atom_charge);
          }
        }

        fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
        if (num_rotamers > 0) {
#pragma omp simd
          for (int i = 0; i < num_nonbonds; ++i) {
            const int &a1 = nonbond_a1_l[i];
            const int &a2 = nonbond_a2_l[i];

            const auto diff_x                = scratch_x_l[a1] - scratch_x_l[a2];
            const auto diff_y                = scratch_y_l[a1] - scratch_y_l[a2];
            const auto diff_z                = scratch_z_l[a1] - scratch_z_l[a2];
            const fp_type distance_two       = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
            const fp_type distance_two_clamp = std::max(distance_two, RMIN_ELEC_SQUARE);
            const fp_type distance           = std::sqrt(distance_two_clamp);

            //  Calculate  Electrostatic  Energy
            const fp_type epsilon =
                mehler_solmajer::A +
                mehler_solmajer::B /
                    (fp_type{1} + mehler_solmajer::rk * std::exp(mehler_solmajer::lambda_B * distance));
            const fp_type r_dielectric = fp_type{1} / (distance * epsilon);
            const fp_type e_elec =
                charge_l[a1] * charge_l[a2] * ELECSCALE * autodock_parameters::coeff_estat * r_dielectric;
            elect_total_eintcal += e_elec;

            // Calcuare desolv
            const fp_type nb_desolv = (vol_l[a2] * (solpar_l[a1] + qsolpar * std::fabs(charge_l[a1])) +
                                       vol_l[a1] * (solpar_l[a2] + qsolpar * std::fabs(charge_l[a2])));

            const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                     std::exp(fp_type{-0.5} / (sigma_square) *distance_two_clamp) * nb_desolv;
            dmap_total_eintcal += e_desolv;
            fp_type e_vdW_Hb{0};
            if (distance_two_clamp < nbc2) {
              //  Find internal energy parameters, i.e.  epsilon and r-equilibrium values...
              //  Lennard-Jones and Hydrogen Bond Potentials
              // This can be precomputed as in intnbtable.cc
              const int xA = xA_default;
              const int xB = nonbond_xB_l[i];

              if (xA != xB) {
                const fp_type cA = nonbond_cA_l[i];
                const fp_type cB = nonbond_cB_l[i];

                const auto log_distance = std::log(distance);
                const fp_type rA        = std::exp(static_cast<fp_type>(xA) * log_distance);
                const fp_type rB        = std::exp(static_cast<fp_type>(xB) * log_distance);

                e_vdW_Hb = std::min(EINTCLAMP, (cA / rA - cB / rB));
              }
            }
            emap_total_eintcal += e_vdW_Hb;
          }
        }
        const fp_type tors_free_energy = static_cast<fp_type>(num_rotamers) * autodock_parameters::coeff_tors;

        const fp_type total_trilinear = emap_total_trilinear + elect_total_trilinear + dmap_total_trilinear;
        const fp_type total_eintcal   = emap_total_eintcal + elect_total_eintcal + dmap_total_eintcal;
        scores_l[scores_index]        = total_trilinear + total_eintcal + tors_free_energy;
      }
    }
    }
  } // namespace

  inline void calc_gradient(const int batch_atoms,
                            const int batch_ligands,
                            const int individuals_per_ligand,
                            const fp_type *__restrict__ x_scratch_b,
                            const fp_type *__restrict__ y_scratch_b,
                            const fp_type *__restrict__ z_scratch_b,
                            const chromosome* __restrict__ chromosomes_b,
                            const int* __restrict__ ligand_fragments_b,
                            const int* __restrict__ ligand_fragments_start_b,
                            const int* __restrict__ frag_indices_start_b,
                            const int* __restrict__ frag_start_indices_b,
                            const int* __restrict__ frag_stop_indices_b,
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
                            gradient *__restrict__ gradients_b,
                            int *__restrict__ active_individuals_b) {
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      const int atom_stride                  = ligand_index * batch_atoms;
      const int num_atoms                    = num_atoms_b[ligand_index];
      const int num_nonbonds                 = num_nonbonds_b[ligand_index + 1] - num_nonbonds_b[ligand_index];
      const int num_rotamers                 = num_rotamers_b[ligand_index];
      const int batch_rotamers               = batch_atoms - 3; // TODO L why is this the max number of rotamers?
      const fp_type *__restrict__ scratch_x  = x_scratch_b + atom_stride * individuals_per_ligand;
      const fp_type *__restrict__ scratch_y  = y_scratch_b + atom_stride * individuals_per_ligand;
      const fp_type *__restrict__ scratch_z  = z_scratch_b + atom_stride * individuals_per_ligand;
      const chromosome* __restrict__ ligand_chromosomes = chromosomes_b + ligand_index * individuals_per_ligand;
      const fp_type *__restrict__ vol_l      = vols_b + atom_stride;
      const fp_type *__restrict__ solpar_l   = solpars_b + atom_stride;
      const fp_type *__restrict__ charge_l   = charges_b + atom_stride;
      const int *__restrict__ map_offsets_l  = map_offsets_b + atom_stride;
      const int *__restrict__ nonbond_a1_l   = nonbond_a1_b + num_nonbonds_b[ligand_index];
      const int *__restrict__ nonbond_a2_l   = nonbond_a2_b + num_nonbonds_b[ligand_index];
      const fp_type *nonbond_cA_l            = nonbond_cA_b + num_nonbonds_b[ligand_index];
      const fp_type *nonbond_cB_l            = nonbond_cB_b + num_nonbonds_b[ligand_index];
      const int *nonbond_xB_l                = nonbond_xB_b + num_nonbonds_b[ligand_index];
      const int* fragments                   = ligand_fragments_b + ligand_fragments_start_b[ligand_index];
      const int* frag_start_indices          = frag_start_indices_b + frag_indices_start_b[ligand_index];
      const int* frag_stop_indices           = frag_stop_indices_b + frag_indices_start_b[ligand_index];
      gradient *__restrict__ gradients_l     = gradients_b + ligand_index * individuals_per_ligand;
      int *__restrict__ active_individuals_l = active_individuals_b + ligand_index * individuals_per_ligand;
      
      std::vector<point3D> dE_dX(batch_atoms);
      std::vector<fp_type> grad(6 + batch_rotamers);         // gradient = dE/dx, dE/dy, dE/dz, dE/dalpha, dE/dbeta, dE/dgamma, dE/d_tors_1, ..., dE/d_tors_n

      for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {
        
        if (!active_individuals_l[individual_index]){
          continue;
        }

        // Reinitialize to 0s all gradient components vectors
        std::fill(dE_dX.begin(), dE_dX.end(), 0);
        std::fill(grad.begin(),  grad.end(),  0);

        const fp_type *__restrict__ scratch_x_l = scratch_x + individual_index * batch_atoms;
        const fp_type *__restrict__ scratch_y_l = scratch_y + individual_index * batch_atoms;
        const fp_type *__restrict__ scratch_z_l = scratch_z + individual_index * batch_atoms;
        
        const fp_type *electro_map = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
        const fp_type *desolv_map  = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

#pragma omp simd
        for (int index = 0; index < num_atoms; ++index) {
          fp_type coord[3]{scratch_x_l[index], scratch_y_l[index], scratch_z_l[index]};

          const auto diff_x = coord[0] - center[0];
          const auto diff_y = coord[1] - center[1];
          const auto diff_z = coord[2] - center[2];

          if (coord[0] < minimum[0] || coord[0] > maximum[0] || coord[1] < minimum[1] ||
              coord[1] > maximum[1] || coord[2] < minimum[2] || coord[2] > maximum[2]) {
            const fp_type penalty_factor = 2 * 2 * ENERGYPENALTY;
            dE_dX[index].x() += penalty_factor * diff_x;
            dE_dX[index].y() += penalty_factor * diff_y;
            dE_dX[index].z() += penalty_factor * diff_z;
          } else {
            const auto &atom_charge = charge_l[index];
            const fp_type *atom_map = grid_maps + map_offsets_l[index];

            coord[0] = (coord[0] - minimum[0]) * inv_spacing;
            coord[1] = (coord[1] - minimum[1]) * inv_spacing;
            coord[2] = (coord[2] - minimum[2]) * inv_spacing;

            const int u0      = static_cast<int>(coord[0]);
            const fp_type p0u = coord[0] - static_cast<fp_type>(u0);
            const fp_type p1u = fp_type{1} - p0u;

            const int v0      = static_cast<int>(coord[1]);
            const fp_type p0v = coord[1] - static_cast<fp_type>(v0);
            const fp_type p1v = fp_type{1} - p0v;

            const int w0      = static_cast<int>(coord[2]);
            const fp_type p0w = coord[2] - static_cast<fp_type>(w0);
            const fp_type p1w = fp_type{1} - p0w;

            // Precompute flattened indices
            const int base_index = FLATTENED_3D(u0, v0, w0, map_index_x, map_index_xy);
            // Trilinear Interpolationp
            fp_type grid_values[8];
            get_grid_values(electro_map + base_index, map_index_x, map_index_xy, grid_values);
            const fp_type constant_factor_el = inv_spacing * atom_charge;
            dE_dX[index].x() += constant_factor_el * (p1w * (p1v * (grid_values[4] - grid_values[0]) + p0v * (grid_values[6] - grid_values[2])) + p0w * (p1v * (grid_values[5] - grid_values[1]) + p0v * (grid_values[7] - grid_values[3])));
            dE_dX[index].y() += constant_factor_el * (p1w * (p1u * (grid_values[2] - grid_values[0]) + p0u * (grid_values[6] - grid_values[4])) + p0w * (p1u * (grid_values[3] - grid_values[1]) + p0u * (grid_values[7] - grid_values[5])));
            dE_dX[index].z() += constant_factor_el * (p1v * (p1u * (grid_values[1] - grid_values[0]) + p0u * (grid_values[5] - grid_values[4])) + p0v * (p1u * (grid_values[3] - grid_values[2]) + p0u * (grid_values[7] - grid_values[6])));

            get_grid_values(atom_map + base_index, map_index_x, map_index_xy, grid_values);
            dE_dX[index].x() += inv_spacing * (p1w * (p1v * (grid_values[4] - grid_values[0]) + p0v * (grid_values[6] - grid_values[2])) + p0w * (p1v * (grid_values[5] - grid_values[1]) + p0v * (grid_values[7] - grid_values[3])));
            dE_dX[index].y() += inv_spacing * (p1w * (p1u * (grid_values[2] - grid_values[0]) + p0u * (grid_values[6] - grid_values[4])) + p0w * (p1u * (grid_values[3] - grid_values[1]) + p0u * (grid_values[7] - grid_values[5])));
            dE_dX[index].z() += inv_spacing * (p1v * (p1u * (grid_values[1] - grid_values[0]) + p0u * (grid_values[5] - grid_values[4])) + p0v * (p1u * (grid_values[3] - grid_values[2]) + p0u * (grid_values[7] - grid_values[6])));            
            
            get_grid_values(desolv_map + base_index, map_index_x, map_index_xy, grid_values);
            const fp_type constant_factor_des = inv_spacing * std::fabs(atom_charge);
            dE_dX[index].x() += constant_factor_des * (p1w * (p1v * (grid_values[4] - grid_values[0]) + p0v * (grid_values[6] - grid_values[2])) + p0w * (p1v * (grid_values[5] - grid_values[1]) + p0v * (grid_values[7] - grid_values[3])));
            dE_dX[index].y() += constant_factor_des * (p1w * (p1u * (grid_values[2] - grid_values[0]) + p0u * (grid_values[6] - grid_values[4])) + p0w * (p1u * (grid_values[3] - grid_values[1]) + p0u * (grid_values[7] - grid_values[5])));
            dE_dX[index].z() += constant_factor_des * (p1v * (p1u * (grid_values[1] - grid_values[0]) + p0u * (grid_values[5] - grid_values[4])) + p0v * (p1u * (grid_values[3] - grid_values[2]) + p0u * (grid_values[7] - grid_values[6])));

            
          }

        }

        if (num_rotamers > 0) {
#pragma omp simd
          for (int i = 0; i < num_nonbonds; ++i) {
            const int &a1 = nonbond_a1_l[i];
            const int &a2 = nonbond_a2_l[i];

            const auto diff_x                = scratch_x_l[a1] - scratch_x_l[a2];
            const auto diff_y                = scratch_y_l[a1] - scratch_y_l[a2];
            const auto diff_z                = scratch_z_l[a1] - scratch_z_l[a2];
            const fp_type distance_two       = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
            const fp_type distance_two_clamp = std::max(distance_two, RMIN_ELEC_SQUARE);
            const fp_type distance           = std::sqrt(distance_two_clamp);

            const fp_type inv_r = fp_type{1} / distance;

            const fp_type dir_x = diff_x * inv_r;
            const fp_type dir_y = diff_y * inv_r;
            const fp_type dir_z = diff_z * inv_r;

            // Electrostatic derivative
            const fp_type exp_term = std::exp(mehler_solmajer::lambda_B * distance);

            const fp_type denom = (fp_type{1} + mehler_solmajer::rk * exp_term);

            const fp_type epsilon = mehler_solmajer::A + mehler_solmajer::B / denom;

            const fp_type epsilon_prime = -mehler_solmajer::B * mehler_solmajer::rk * mehler_solmajer::lambda_B * exp_term / (denom * denom);

            const fp_type C = charge_l[a1] * charge_l[a2] * ELECSCALE * autodock_parameters::coeff_estat;

            const fp_type f = distance * epsilon;
            const fp_type f_prime = epsilon + distance * epsilon_prime;

            const fp_type dE_dr_elec = -C * f_prime / (f * f);

            // Calcuare desolv
            const fp_type nb_desolv = (vol_l[a2] * (solpar_l[a1] + qsolpar * std::fabs(charge_l[a1])) +
                                       vol_l[a1] * (solpar_l[a2] + qsolpar * std::fabs(charge_l[a2])));

            const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                     std::exp(fp_type{-0.5} / (sigma_square) *distance_two_clamp) * nb_desolv;

            // Desolvation derivative
            const fp_type dE_dr_desolv = (-distance / sigma_square) * e_desolv;


            fp_type dE_dr_vdw = 0;
            if (distance_two_clamp < nbc2) {
              const int xA = xA_default;
              const int xB = nonbond_xB_l[i];

              if (xA != xB) {
                const fp_type cA = nonbond_cA_l[i];
                const fp_type cB = nonbond_cB_l[i];

                // r^{-x}
                const auto log_distance = std::log(distance);
                const fp_type rA        = std::exp(-static_cast<fp_type>(xA) * log_distance);
                const fp_type rB        = std::exp(-static_cast<fp_type>(xB) * log_distance);
                
                // VdW derivative
                const fp_type rA1 = rA * inv_r; // r^{-(xA+1)}
                const fp_type rB1 = rB * inv_r; // r^{-(xB+1)}
                dE_dr_vdw = -static_cast<fp_type>(xA) * cA * rA1 + static_cast<fp_type>(xB) * cB * rB1;
              }
            }

            const fp_type dE_dr = dE_dr_elec + dE_dr_desolv + dE_dr_vdw;

            // Accumulate into dE_dX
            const fp_type gx = dE_dr * dir_x;
            const fp_type gy = dE_dr * dir_y;
            const fp_type gz = dE_dr * dir_z;

            dE_dX[a1].x() += gx;
            dE_dX[a1].y() += gy;
            dE_dX[a1].z() += gz;

            dE_dX[a2].x() -= gx;
            dE_dX[a2].y() -= gy;
            dE_dX[a2].z() -= gz;

          }
        }

        point3D ligand_COM{fp_type{0}};
        for (int i = 0; i < num_atoms; ++i) {
          ligand_COM.x() += scratch_x_l[i];
          ligand_COM.y() += scratch_y_l[i];
          ligand_COM.z() += scratch_z_l[i];
        }
        ligand_COM.x() /= static_cast<fp_type>(num_atoms);
        ligand_COM.y() /= static_cast<fp_type>(num_atoms);
        ligand_COM.z() /= static_cast<fp_type>(num_atoms);

        point3D tau{fp_type{0}};
        for (int i = 0; i < num_atoms; ++i) {
          point3D r;
          r.x() = scratch_x_l[i] - ligand_COM.x();
          r.y() = scratch_y_l[i] - ligand_COM.y();
          r.z() = scratch_z_l[i] - ligand_COM.z();
          const point3D torque = r.cross(dE_dX[i]);
          tau.x() += torque.x();
          tau.y() += torque.y();
          tau.z() += torque.z();
        }

        const chromosome &chrom = ligand_chromosomes[individual_index];
        const fp_type rad_beta = deg_to_rad(chrom[4]);
        const fp_type rad_gamma = deg_to_rad(chrom[5]);

        const fp_type sin_beta = std::sin(rad_beta);
        const fp_type cos_beta = std::cos(rad_beta);
        const fp_type sin_gamma = std::sin(rad_gamma);
        const fp_type cos_gamma = std::cos(rad_gamma);

        point3D u_gamma;
        u_gamma.x() = fp_type{0};
        u_gamma.y() = fp_type{0};
        u_gamma.z() = fp_type{1};

        point3D u_beta;
        u_beta.x() = -sin_gamma;
        u_beta.y() = cos_gamma;
        u_beta.z() = fp_type{0};

        point3D u_alpha;
        u_alpha.x() = cos_beta * cos_gamma;
        u_alpha.y() = cos_beta * sin_gamma;
        u_alpha.z() = -sin_beta;

        grad[3] = tau.inner_product(u_alpha);
        grad[4] = tau.inner_product(u_beta);
        grad[5] = tau.inner_product(u_gamma);

        for (int index = 0; index < num_atoms; ++index) {
          grad[0] += dE_dX[index].x();
          grad[1] += dE_dX[index].y();
          grad[2] += dE_dX[index].z();
        }

        // Torsion derivatives
        for (int t = 0; t < num_rotamers; ++t) {
          const int* frag_mask = fragments + t * num_atoms;

          const int a1 = frag_start_indices[t];
          const int a2 = frag_stop_indices[t];

          point3D axis;
          axis.x() = scratch_x_l[a2] - scratch_x_l[a1];
          axis.y() = scratch_y_l[a2] - scratch_y_l[a1];
          axis.z() = scratch_z_l[a2] - scratch_z_l[a1];

          const fp_type norm = std::sqrt(axis.x() * axis.x() + axis.y() * axis.y() + axis.z() * axis.z());

          const fp_type inv_norm = fp_type{1} / norm;
          axis.x() *= inv_norm;
          axis.y() *= inv_norm;
          axis.z() *= inv_norm;

          point3D r;
#pragma omp simd
          for (int i = 0; i < num_atoms; ++i) {
            if (frag_mask[i] != 0) {
              r.x() = scratch_x_l[i] - scratch_x_l[a1];
              r.y() = scratch_y_l[i] - scratch_y_l[a1];
              r.z() = scratch_z_l[i] - scratch_z_l[a1];
              grad[6 + t] += dE_dX[i].inner_product(axis.cross(r));
            }
          }
        }

        // Write gradient on the buffer
        gradient &grad_l = gradients_l[individual_index];
        for (int i = 0; i < 6 + num_rotamers; ++i) {
          grad_l[i] = grad[i];
        }
        
      }
    }
  };


  template<>
  void adt_score_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<adt_region_name>(calc_energy,
                                      batch_atoms,
                                      batch_ligands,
                                      scores_per_ligand,
                                      x_scratch_b,
                                      y_scratch_b,
                                      z_scratch_b,
                                      vols_b,
                                      solpars_b,
                                      charges_b,
                                      num_atoms_b,
                                      num_rotamers_b,
                                      num_nonbonds_b,
                                      nonbond_a1_b,
                                      nonbond_a2_b,
                                      nonbond_cA_b,
                                      nonbond_cB_b,
                                      nonbond_xB_b,
                                      grid_maps,
                                      minimum,
                                      maximum,
                                      center,
                                      map_offsets_b,
                                      map_index_x,
                                      map_index_xy,
                                      map_index_xyz,
                                      scores_b);
  }

  template<>
  void adt_gradient_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<this->adt_region_name>(calc_gradient,
                                            batch_atoms,
                                            batch_ligands,
                                            scores_per_ligand,
                                            x_scratch_b,
                                            y_scratch_b,
                                            z_scratch_b,
                                            chromosomes_b,
                                            ligand_fragments_b,
                                            ligand_fragments_start_b,
                                            frag_indices_start_b,
                                            frag_start_indices_b,
                                            frag_stop_indices_b,
                                            vols_b,
                                            solpars_b,
                                            charges_b,
                                            num_atoms_b,
                                            num_rotamers_b,
                                            num_nonbonds_b,
                                            nonbond_a1_b,
                                            nonbond_a2_b,
                                            nonbond_cA_b,
                                            nonbond_cB_b,
                                            nonbond_xB_b,
                                            grid_maps,
                                            minimum,
                                            maximum,
                                            center,
                                            map_offsets_b,
                                            map_index_x,
                                            map_index_xy,
                                            map_index_xyz,
                                            gradients_b,
                                            active_individuals_b
    );
  }
} // namespace mudock