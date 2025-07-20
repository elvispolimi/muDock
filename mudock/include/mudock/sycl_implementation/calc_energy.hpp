#pragma once

#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/type_alias.hpp>
#include <sycl/sycl.hpp>

#define FLATTENED_3D(x, y, z, index_x, index_xy) (index_xy * (z) + (y) * index_x + (x))

namespace mudock {
  inline fp_type trilinear_interpolation_sycl(const fp_type* __restrict__ map,
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

  template<int MAX_ATOMS>
  inline fp_type calc_intra_energy(const fp_type* ligand_x,
                                   const fp_type* ligand_y,
                                   const fp_type* ligand_z,
                                   const fp_type* ligand_vol,
                                   const fp_type* ligand_solpar,
                                   const fp_type* ligand_charge,
                                   const int num_atoms,
                                   const int num_rotamers,
                                   const int ligand_num_nonbonds,
                                   const int* __restrict__ ligand_nonbond_a1,
                                   const int* __restrict__ ligand_nonbond_a2,
                                   const fp_type* __restrict__ ligand_nonbond_cA,
                                   const fp_type* __restrict__ ligand_nonbond_cB,
                                   const int* __restrict__ ligand_nonbond_xB,
                                   const point3D minimum,
                                   const point3D maximum,
                                   const point3D center,
                                   const int map_index_x,
                                   const int map_index_xy,
                                   const int map_index_xyz,
                                   const fp_type* __restrict__ grid_maps,
                                   const int* __restrict__ map_ligand_offsets,
                                   sycl::nd_item<1> it) {
    const int workgroup_size       = it.get_local_range(0);
    const int workgroup_id         = it.get_group(0);
    const int workitem_id_in_group = it.get_local_id(0);

    const fp_type* electro_map = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
    const fp_type* desolv_map  = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

    // Calculate energy
    fp_type elect_total_trilinear = 0;
    fp_type emap_total_trilinear  = 0;
    fp_type dmap_total_trilinear  = 0;
    for (int atom_index = workitem_id_in_group; atom_index < MAX_ATOMS; atom_index += workgroup_size) {
      if (atom_index < num_atoms) {
        fp_type coord_tex[3]{ligand_x[atom_index], ligand_y[atom_index], ligand_z[atom_index]};
        if (coord_tex[0] < minimum.x() || coord_tex[0] > maximum.x() || coord_tex[1] < minimum.y() ||
            coord_tex[1] > maximum.y() || coord_tex[2] < minimum.z() || coord_tex[2] > maximum.z()) {
          // Is outside
          const auto diff_x          = coord_tex[0] - center.x();
          const auto diff_y          = coord_tex[1] - center.y();
          const auto diff_z          = coord_tex[2] - center.z();
          const fp_type distance_two = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;

          const fp_type epenalty = distance_two * ENERGYPENALTY;
          elect_total_trilinear += epenalty;
          emap_total_trilinear += epenalty;
        } else {
          // Is inside
          // Center atom coordinates on the grid center
          coord_tex[0]            = (coord_tex[0] - minimum.x()) * inv_spacing,
          coord_tex[1]            = (coord_tex[1] - minimum.y()) * inv_spacing;
          coord_tex[2]            = (coord_tex[2] - minimum.z()) * inv_spacing;
          const auto& charge      = ligand_charge[atom_index];
          const fp_type* atom_map = grid_maps + map_ligand_offsets[atom_index];

          const int u0      = coord_tex[0];
          const fp_type p0u = coord_tex[0] - static_cast<fp_type>(u0);
          const fp_type p1u = fp_type{1} - p0u;

          const int v0      = coord_tex[1];
          const fp_type p0v = coord_tex[1] - static_cast<fp_type>(v0);
          const fp_type p1v = fp_type{1} - p0v;

          const int w0      = coord_tex[2];
          const fp_type p0w = coord_tex[2] - static_cast<fp_type>(w0);
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
              trilinear_interpolation_sycl(electro_map + base_index, coeffs, map_index_x, map_index_xy) *
              charge;
          emap_total_trilinear +=
              trilinear_interpolation_sycl(atom_map + base_index, coeffs, map_index_x, map_index_xy);
          dmap_total_trilinear +=
              trilinear_interpolation_sycl(desolv_map + base_index, coeffs, map_index_x, map_index_xy) *
              std::fabs(charge);
        }
      }
    }

    fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
    if (num_rotamers > 0)
      for (int nonbond_list = it.get_local_id(0); nonbond_list < ligand_num_nonbonds;
           nonbond_list += it.get_local_range(0)) {
        const int& a1 = ligand_nonbond_a1[nonbond_list];
        const int& a2 = ligand_nonbond_a2[nonbond_list];

        const auto diff_x                = ligand_x[a1] - ligand_x[a2];
        const auto diff_y                = ligand_y[a1] - ligand_y[a2];
        const auto diff_z                = ligand_z[a1] - ligand_z[a2];
        const fp_type distance_two       = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
        const fp_type distance_two_clamp = sycl::max(distance_two, RMIN_ELEC_SQUARE);
        const fp_type distance           = sycl::sqrt(distance_two_clamp);

        //  Calculate  Electrostatic  Energy
        const fp_type epsilon =
            mehler_solmajer::A +
            mehler_solmajer::B /
                (fp_type{1} + mehler_solmajer::rk * sycl::exp(mehler_solmajer::lambda_B * distance));
        const fp_type r_dielectric = fp_type{1} / (distance * epsilon);
        const fp_type e_elec       = ligand_charge[a1] * ligand_charge[a2] * ELECSCALE *
                               autodock_parameters::coeff_estat * r_dielectric;
        elect_total_eintcal += e_elec;

        // Calcuate desolv
        const fp_type nb_desolv =
            (ligand_vol[a2] * (ligand_solpar[a1] + qsolpar * sycl::fabs(ligand_charge[a1])) +
             ligand_vol[a1] * (ligand_solpar[a2] + qsolpar * sycl::fabs(ligand_charge[a2])));

        const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                 sycl::exp(fp_type{-0.5} / (sigma * sigma) * distance_two_clamp) * nb_desolv;
        dmap_total_eintcal += e_desolv;

        fp_type e_vdW_Hb{0};
        if (distance_two_clamp < nbc2) {
          //  Find internal energy parameters, i.e.  epsilon and r-equilibrium values...
          //  Lennard-Jones and Hydrogen Bond Potentials
          // This can be precomputed as in intnbtable.cc
          const int xA = xA_default;
          const int xB = ligand_nonbond_xB[nonbond_list];

          if (xA != xB) {
            const fp_type cA = ligand_nonbond_cA[nonbond_list];
            const fp_type cB = ligand_nonbond_cB[nonbond_list];

            const auto log_distance = sycl::log(distance);
            const fp_type rA        = sycl::exp(static_cast<fp_type>(xA) * log_distance);
            const fp_type rB        = sycl::exp(static_cast<fp_type>(xB) * log_distance);

            e_vdW_Hb = sycl::min(EINTCLAMP, (cA / rA - cB / rB));
          }
        }
        emap_total_eintcal += e_vdW_Hb;
      }
    return elect_total_eintcal + emap_total_eintcal + dmap_total_eintcal + elect_total_trilinear +
           dmap_total_trilinear + emap_total_trilinear;
  };
} // namespace mudock
