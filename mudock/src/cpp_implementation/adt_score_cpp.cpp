#include <cstring>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/cpp_implementation/adt_score_cpp.hpp>

#define FLATTENED_2D(x, y, index_x)              ((y) * index_x + (x))
#define FLATTENED_3D(x, y, z, index_x, index_xy) (index_xy * (z) + (y) * index_x + (x))

namespace mudock {
  inline fp_type trilinear_interpolation(const fp_type *__restrict__ map,
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
                              const fp_type* (&out_values)[8]) {
    out_values[0] = &map[0];
    out_values[1] = &map[map_index_xy];
    out_values[2] = &map[map_index_x];
    out_values[3] = &map[map_index_x + map_index_xy];
    out_values[4] = &map[1];
    out_values[5] = &map[1 + map_index_xy];
    out_values[6] = &map[1 + map_index_x];
    out_values[7] = &map[1 + map_index_x + map_index_xy];
  }

  inline fp_type dot_product(const point3D& u, const point3D& v){
    fp_type result;
    result = u.x() * v.x() + u.y() * v.y() + u.z() * v.z();
    return result;
  }

  inline void calc_energy(const int batch_atoms,
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

            const int u0      = coord[0];
            const fp_type p0u = coord[0] - static_cast<fp_type>(u0);
            const fp_type p1u = fp_type{1} - p0u;

            const int v0      = coord[1];
            const fp_type p0v = coord[1] - static_cast<fp_type>(v0);
            const fp_type p1v = fp_type{1} - p0v;

            const int w0      = coord[2];
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
        const fp_type tors_free_energy = num_rotamers * autodock_parameters::coeff_tors;

        const fp_type total_trilinear = emap_total_trilinear + elect_total_trilinear + dmap_total_trilinear;
        const fp_type total_eintcal   = emap_total_eintcal + elect_total_eintcal + dmap_total_eintcal;
        scores_l[scores_index]        = total_trilinear + total_eintcal + tors_free_energy;
      }
    }
  };

  inline void calc_gradient(const int batch_atoms,
                          const int batch_ligands,
                          const int individuals_per_ligand,
                          const fp_type *__restrict__ x_scratch_b,
                          const fp_type *__restrict__ y_scratch_b,
                          const fp_type *__restrict__ z_scratch_b,
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
                            fp_type *__restrict__ gradients_b) {
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

      const int* fragments = ligand_fragments_b + ligand_fragments_start_b[ligand_index];
      const int* frag_start_indices = frag_start_indices_b + frag_indices_start_b[ligand_index];
      const int* frag_stop_indices = frag_stop_indices_b + frag_indices_start_b[ligand_index];

      fp_type *__restrict__ gradients_l = gradients_b + ligand_index * individuals_per_ligand;
      for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {
        const fp_type *__restrict__ scratch_x_l = scratch_x + individual_index * batch_atoms;
        const fp_type *__restrict__ scratch_y_l = scratch_y + individual_index * batch_atoms;
        const fp_type *__restrict__ scratch_z_l = scratch_z + individual_index * batch_atoms;
        
        std::vector<point3D> dE_dX(num_atoms);
        std::vector<point3D> dX_dalpha(num_atoms);
        std::vector<point3D> dX_dbeta(num_atoms);
        std::vector<point3D> dX_dgamma(num_atoms);
        std::vector<std::vector<point3D>> dX_drot(num_rotamers, std::vector<point3D>(num_atoms));
        
        // gradient = dE/dx, dE/dy, dE/dz, dE/dalpha, dE/dbeta, dE/dgamma, dE/d_tors_1, ..., dE/d_tors_n
        std::vector<fp_type> grad(6 + num_rotamers, 0);

        const fp_type *electro_map = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
        const fp_type *desolv_map  = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

#pragma omp simd
        for (int index = 0; index < num_atoms; ++index) {
          fp_type coord[3]{scratch_x_l[index], scratch_y_l[index], scratch_z_l[index]};

          const auto diff_x      = coord[0] - center[0];
          const auto diff_y      = coord[1] - center[1];
          const auto diff_z      = coord[2] - center[2];

          if (coord[0] < minimum[0] || coord[0] > maximum[0] || coord[1] < minimum[1] ||
              coord[1] > maximum[1] || coord[2] < minimum[2] || coord[2] > maximum[2]) {
            const fp_type dist     = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
            const fp_type epenalty = dist * ENERGYPENALTY;
            // TODO gestire calcolo gradiente in questo if
            // dE_dX.x() +=
            // dE_dX.y() +=
            // dE_dX.z() +=
          } else {
            const auto &atom_charge = charge_l[index];
            const fp_type *atom_map = grid_maps + map_offsets_l[index];

            coord[0] = (coord[0] - minimum[0]) * inv_spacing;
            coord[1] = (coord[1] - minimum[1]) * inv_spacing;
            coord[2] = (coord[2] - minimum[2]) * inv_spacing;

            const int u0      = coord[0];
            const fp_type p0u = coord[0] - static_cast<fp_type>(u0);
            const fp_type p1u = fp_type{1} - p0u;

            const int v0      = coord[1];
            const fp_type p0v = coord[1] - static_cast<fp_type>(v0);
            const fp_type p1v = fp_type{1} - p0v;

            const int w0      = coord[2];
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
            const fp_type* grid_values[8];
            get_grid_values(electro_map + base_index, map_index_x, map_index_xy, grid_values);
            // TODO non penso ci sia bisogno di moltiplicare per inv_spacing perché viene già fatto quando normalizza coord[]? controllare. LA QUESTIONE VALE PER TUTTE E TRE LE ENERGIE
            dE_dX[index].x() += atom_charge * (p1w * (p1v * (grid_values[4] - grid_values[0]) + p0v * (grid_values[6] - grid_values[2])) + p0w * (p1v * (grid_values[5] - grid_values[1]) + p0v * (grid_values[7] - grid_values[3])));
            dE_dX[index].y() += atom_charge * (p1w * (p1u * (grid_values[2] - grid_values[0]) + p0u * (grid_values[6] - grid_values[4])) + p0w * (p1u * (grid_values[3] - grid_values[1]) + p0u * (grid_values[7] - grid_values[5])));
            dE_dX[index].z() += atom_charge * (p1v * (p1u * (grid_values[1] - grid_values[0]) + p0u * (grid_values[5] - grid_values[4])) + p0v * (p1u * (grid_values[3] - grid_values[2]) + p0u * (grid_values[7] - grid_values[6])));

            get_grid_values(atom_map + base_index, map_index_x, map_index_xy, grid_values);
            dE_dX[index].x() += (p1w * (p1v * (grid_values[4] - grid_values[0]) + p0v * (grid_values[6] - grid_values[2])) + p0w * (p1v * (grid_values[5] - grid_values[1]) + p0v * (grid_values[7] - grid_values[3])));
            dE_dX[index].y() += (p1w * (p1u * (grid_values[2] - grid_values[0]) + p0u * (grid_values[6] - grid_values[4])) + p0w * (p1u * (grid_values[3] - grid_values[1]) + p0u * (grid_values[7] - grid_values[5])));
            dE_dX[index].z() += (p1v * (p1u * (grid_values[1] - grid_values[0]) + p0u * (grid_values[5] - grid_values[4])) + p0v * (p1u * (grid_values[3] - grid_values[2]) + p0u * (grid_values[7] - grid_values[6])));            
            
            get_grid_values(desolv_map + base_index, map_index_x, map_index_xy, grid_values);
            dE_dX[index].x() += std::fabs(atom_charge) * (p1w * (p1v * (grid_values[4] - grid_values[0]) + p0v * (grid_values[6] - grid_values[2])) + p0w * (p1v * (grid_values[5] - grid_values[1]) + p0v * (grid_values[7] - grid_values[3])));
            dE_dX[index].y() += std::fabs(atom_charge) * (p1w * (p1u * (grid_values[2] - grid_values[0]) + p0u * (grid_values[6] - grid_values[4])) + p0w * (p1u * (grid_values[3] - grid_values[1]) + p0u * (grid_values[7] - grid_values[5])));
            dE_dX[index].z() += std::fabs(atom_charge) * (p1v * (p1u * (grid_values[1] - grid_values[0]) + p0u * (grid_values[5] - grid_values[4])) + p0v * (p1u * (grid_values[3] - grid_values[2]) + p0u * (grid_values[7] - grid_values[6])));

            // Compute dX/d_rot
            // alpha -> rotation on z-axis
            // TODO check if it is correct to use diff_x/y/z
            // TODO check if it is fine to consider infinitesimal rotations so we can just use the versors as u or we shuold consider the composition of rotations in order
            const point3D x_i = point3D{diff_x, diff_y, diff_z};

            const point3D z_axis = point3D{fp_type{0}, fp_type{0}, fp_type{1}};
            const point3D y_axis = point3D{fp_type{0}, fp_type{1}, fp_type{0}};
            const point3D x_axis = point3D{fp_type{1}, fp_type{0}, fp_type{0}};

            dX_dalpha[index] = z_axis.cross(x_i);
            dX_dbeta[index]  = y_axis.cross(x_i);
            dX_dgamma[index] = x_axis.cross(x_i);

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

            const fp_type inv_r = fp_type{1} / distance;

            const fp_type dir_x = diff_x * inv_r;
            const fp_type dir_y = diff_y * inv_r;
            const fp_type dir_z = diff_z * inv_r;

            // Electrostatic derivative
            const fp_type exp_term =
                std::exp(mehler_solmajer::lambda_B * distance);

            const fp_type denom =
                (fp_type{1} + mehler_solmajer::rk * exp_term);

            const fp_type epsilon =
                mehler_solmajer::A +
                mehler_solmajer::B / denom;

            const fp_type epsilon_prime =
                -mehler_solmajer::B *
                mehler_solmajer::rk *
                mehler_solmajer::lambda_B *
                exp_term / (denom * denom);

            const fp_type C =
                charge_l[a1] * charge_l[a2] *
                ELECSCALE * autodock_parameters::coeff_estat;

            const fp_type f = distance * epsilon;
            const fp_type f_prime = epsilon + distance * epsilon_prime;

            const fp_type dE_dr_elec = -C * f_prime / (f * f);

            // Calcuare desolv
            const fp_type nb_desolv = (vol_l[a2] * (solpar_l[a1] + qsolpar * std::fabs(charge_l[a1])) +
                                       vol_l[a1] * (solpar_l[a2] + qsolpar * std::fabs(charge_l[a2])));

            const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                     std::exp(fp_type{-0.5} / (sigma_square) *distance_two_clamp) * nb_desolv;
            dmap_total_eintcal += e_desolv;

            // Desolvation derivative
            const fp_type dE_dr_desolv = (-distance / sigma_square) * e_desolv;


            fp_type dE_dr_vdw = 0;
            if (distance_two_clamp < nbc2) {
              //  Find internal energy parameters, i.e.  epsilon and r-equilibrium values...
              //  Lennard-Jones and Hydrogen Bond Potentials
              // This can be precomputed as in intnbtable.cc
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
                dE_dr_vdw = -xA * cA * rA1 + xB * cB * rB1;
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

          const fp_type inv_norm = 1.0 / norm;
          axis.x() *= inv_norm;
          axis.y() *= inv_norm;
          axis.z() *= inv_norm;

#pragma omp simd
          for (int i = 0; i < num_atoms; ++i) {
            if (frag_mask[i] != 0) {

              point3D r;
              r.x() = scratch_x_l[i] - scratch_x_l[a1];
              r.y() = scratch_y_l[i] - scratch_y_l[a1];
              r.z() = scratch_z_l[i] - scratch_z_l[a1];

              dX_drot[t][i] = axis.cross(r);
            }
          }
        }

        // Accumulate gradient over atoms: dE/dtheta = SUM_i(dE/dX_i * dX_i/dtheta)
        for (int index = 0; index < num_atoms; ++index){
          grad[0] += dE_dX[index].x();
          grad[1] += dE_dX[index].y();
          grad[2] += dE_dX[index].z();

          grad[3] += dot_product(dE_dX[index], dX_dalpha[index]);
          grad[4] += dot_product(dE_dX[index], dX_dbeta[index]);
          grad[5] += dot_product(dE_dX[index], dX_dgamma[index]);
        }

        for (int t = 0; t < num_rotamers; ++t) {
          for (int index = 0; index < num_atoms; ++index) {
            grad[6 + t] += dot_product(dE_dX[index], dX_drot[t][index]);
          }
        }

        for (int i = 0; i < 6 + num_rotamers; ++i) {
          gradients_l[individual_index * (6 + num_rotamers) + i] = grad[i];
        }
        
      }
    }
  };


  template<>
  void adt_score_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<this->adt_region_name>(calc_energy,
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
                                            gradients_b
    );
  }
} // namespace mudock