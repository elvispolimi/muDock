#include <cstring>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
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
} // namespace mudock
