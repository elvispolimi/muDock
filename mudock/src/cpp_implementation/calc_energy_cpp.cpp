#include <cstdio>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/cpp_implementation/calc_energy_cpp.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/type_alias.hpp>

#define FLATTENED_2D(x, y, index_x)              ((y) * index_x + (x))
#define FLATTENED_3D(x, y, z, index_x, index_xy) (index_xy * (z) + (y) * index_x + (x))

namespace mudock {
  inline fp_type trilinear_interpolation(const fp_type* __restrict__ map,
                                         const fp_type* __restrict__ coeffs,
                                         const int& map_index_x,
                                         const int& map_index_xy) {
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

  template<>
  fp_type calc_energy<cpu_vectorization::AUTO>(const fp_type* __restrict__ ligand_x,
                                               const fp_type* __restrict__ ligand_y,
                                               const fp_type* __restrict__ ligand_z,
                                               const fp_type* __restrict__ ligand_vol,
                                               const fp_type* __restrict__ ligand_solpar,
                                               const fp_type* __restrict__ ligand_charge,
                                               const int* __restrict__ map_ligand_offsets,
                                               const int num_atoms,
                                               const int n_torsions,
                                               const int num_nonbond,
                                               const int* __restrict__ non_bond_list_a1,
                                               const int* __restrict__ non_bond_list_a2,
                                               const fp_type* __restrict__ cA_list,
                                               const fp_type* __restrict__ cB_list,
                                               const int* __restrict__ xB_list,
                                               const fp_type* __restrict__ minimum,
                                               const fp_type* __restrict__ maximum,
                                               const fp_type* __restrict__ center,
                                               const int map_index_x,
                                               const int map_index_xy,
                                               const fp_type* __restrict__ grid_maps,
                                               const fp_type* __restrict__ electro_map,
                                               const fp_type* __restrict__ desolv_map) {
    fp_type elect_total_trilinear = 0;
    fp_type emap_total_trilinear  = 0;
    fp_type dmap_total_trilinear  = 0;

#pragma omp simd
    for (int index = 0; index < num_atoms; ++index) {
      fp_type coord[3]{ligand_x[index], ligand_y[index], ligand_z[index]};

      if (coord[0] < minimum[0] || coord[0] > maximum[0] || coord[1] < minimum[1] || coord[1] > maximum[1] ||
          coord[2] < minimum[2] || coord[2] > maximum[2]) {
        const auto diff_x      = coord[0] - center[0];
        const auto diff_y      = coord[1] - center[1];
        const auto diff_z      = coord[2] - center[2];
        const fp_type dist     = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
        const fp_type epenalty = dist * ENERGYPENALTY;
        elect_total_trilinear += epenalty;
        emap_total_trilinear += epenalty;
      } else {
        const auto& atom_charge = ligand_charge[index];
        const fp_type* atom_map = grid_maps + map_ligand_offsets[index];

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
    if (n_torsions > 0) {
#pragma omp simd
      for (int i = 0; i < num_nonbond; ++i) {
        const int& a1 = non_bond_list_a1[i];
        const int& a2 = non_bond_list_a2[i];

        const auto diff_x                = ligand_x[a1] - ligand_x[a2];
        const auto diff_y                = ligand_y[a1] - ligand_y[a2];
        const auto diff_z                = ligand_z[a1] - ligand_z[a2];
        const fp_type distance_two       = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
        const fp_type distance_two_clamp = std::max(distance_two, RMIN_ELEC_SQUARE);
        const fp_type distance           = std::sqrt(distance_two_clamp);

        //  Calculate  Electrostatic  Energy
        const fp_type epsilon =
            mehler_solmajer::A +
            mehler_solmajer::B /
                (fp_type{1} + mehler_solmajer::rk * std::exp(mehler_solmajer::lambda_B * distance));
        const fp_type r_dielectric = fp_type{1} / (distance * epsilon);
        const fp_type e_elec       = ligand_charge[a1] * ligand_charge[a2] * ELECSCALE *
                               autodock_parameters::coeff_estat * r_dielectric;
        elect_total_eintcal += e_elec;

        // Calcuare desolv
        const fp_type nb_desolv =
            (ligand_vol[a2] * (ligand_solpar[a1] + qsolpar * std::fabs(ligand_charge[a1])) +
             ligand_vol[a1] * (ligand_solpar[a2] + qsolpar * std::fabs(ligand_charge[a2])));

        const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                 std::exp(fp_type{-0.5} / (sigma_square) *distance_two_clamp) * nb_desolv;
        dmap_total_eintcal += e_desolv;
        fp_type e_vdW_Hb{0};
        if (distance_two_clamp < nbc2) {
          //  Find internal energy parameters, i.e.  epsilon and r-equilibrium values...
          //  Lennard-Jones and Hydrogen Bond Potentials
          // This can be precomputed as in intnbtable.cc
          const int xA = xA_default;
          const int xB = xB_list[i];

          if (xA != xB) {
            const fp_type cA = cA_list[i];
            const fp_type cB = cB_list[i];

            const auto log_distance = std::log(distance);
            const fp_type rA        = std::exp(static_cast<fp_type>(xA) * log_distance);
            const fp_type rB        = std::exp(static_cast<fp_type>(xB) * log_distance);

            e_vdW_Hb = std::min(EINTCLAMP, (cA / rA - cB / rB));
          }
        }
        emap_total_eintcal += e_vdW_Hb;
      }
    }
    const fp_type tors_free_energy = n_torsions * autodock_parameters::coeff_tors;

    const fp_type total_trilinear = emap_total_trilinear + elect_total_trilinear + dmap_total_trilinear;
    const fp_type total_eintcal   = emap_total_eintcal + elect_total_eintcal + dmap_total_eintcal;
    return total_trilinear + total_eintcal + tors_free_energy;
  }

} // namespace mudock
