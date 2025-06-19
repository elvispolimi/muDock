#pragma once

#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  __device__ static constexpr fp_type EINTCLAMP_CUDA{EINTCLAMP};
  __device__ static constexpr fp_type lambda{0.003627};
  __device__ static constexpr fp_type epsilon0{78.4};
  __device__ static constexpr fp_type A{-8.5525};
  __device__ static constexpr fp_type B = epsilon0 - A;
  __device__ static constexpr fp_type rk{7.7839};
  __device__ static constexpr fp_type lambda_B = -lambda * B;

  __device__ static constexpr fp_type RMIN_ELEC_SQUARE_CUDA = RMIN_ELEC_SQUARE;
  __device__ static constexpr fp_type sigma_square_cuda     = sigma_square;

  extern __device__ __constant__ fp_type map_min_const[3];
  extern __device__ __constant__ fp_type map_max_const[3];
  extern __device__ __constant__ fp_type map_center_const[3];

  __device__ inline fp_type trilinear_interpolation_cuda(const int coord[],
                                                         const cudaTextureObject_t& tex,
                                                         const fp_type* __restrict__ coeffs) {
    // Interpolation CUDA
    fp_type value{0};

    value = coeffs[0] * tex3D<fp_type>(tex, coord[0], coord[1], coord[2]) + value;
    value = coeffs[1] * tex3D<fp_type>(tex, coord[0], coord[1], coord[2] + 1) + value;
    value = coeffs[2] * tex3D<fp_type>(tex, coord[0], coord[1] + 1, coord[2]) + value;
    value = coeffs[3] * tex3D<fp_type>(tex, coord[0], coord[1] + 1, coord[2] + 1) + value;
    value = coeffs[4] * tex3D<fp_type>(tex, coord[0] + 1, coord[1], coord[2]) + value;
    value = coeffs[5] * tex3D<fp_type>(tex, coord[0] + 1, coord[1], coord[2] + 1) + value;
    value = coeffs[6] * tex3D<fp_type>(tex, coord[0] + 1, coord[1] + 1, coord[2]) + value;
    value = coeffs[7] * tex3D<fp_type>(tex, coord[0] + 1, coord[1] + 1, coord[2] + 1) + value;

    return value;
  }

  template<int MAX_ATOMS>
  __device__ inline fp_type calc_intra_energy(const fp_type* ligand_x,
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
                                              const cudaTextureObject_t* __restrict__ atom_textures,
                                              const int* __restrict__ atom_tex_indexes,
                                              const cudaTextureObject_t electro_texture,
                                              const cudaTextureObject_t desolv_texture) {
    // Calculate energy
    fp_type elect_total_trilinear = 0;
    fp_type emap_total_trilinear  = 0;
    fp_type dmap_total_trilinear  = 0;
#pragma unroll
    for (int atom_index = threadIdx.x; atom_index < MAX_ATOMS; atom_index += blockDim.x) {
      if (atom_index < num_atoms) {
        fp_type coord_tex[3]{ligand_x[atom_index], ligand_y[atom_index], ligand_z[atom_index]};

        if (coord_tex[0] < map_min_const[0] || coord_tex[0] > map_max_const[0] ||
            coord_tex[1] < map_min_const[1] || coord_tex[1] > map_max_const[1] ||
            coord_tex[2] < map_min_const[2] || coord_tex[2] > map_max_const[2]) {
          // Is outside
          const auto diff_x          = coord_tex[0] - map_center_const[0];
          const auto diff_y          = coord_tex[1] - map_center_const[1];
          const auto diff_z          = coord_tex[2] - map_center_const[2];
          const fp_type distance_two = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;

          const fp_type epenalty = distance_two * ENERGYPENALTY;
          elect_total_trilinear += epenalty;
          emap_total_trilinear += epenalty;
        } else {
          // Is inside
          // Center atom coordinates on the grid center
          coord_tex[0] = (coord_tex[0] - map_min_const[0]) * inv_spacing,
          coord_tex[1] = (coord_tex[1] - map_min_const[1]) * inv_spacing;
          coord_tex[2] = (coord_tex[2] - map_min_const[2]) * inv_spacing;

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
          const int int_coord[3]  = {u0, v0, w0};
          //  TODO check approximations with in hardware interpolation
          // elect_total_trilinear +=
          //     tex3D<fp_type>(electro_texture, coord_tex[0], coord_tex[1], coord_tex[2]) *
          //     l_ligand_charge[atom_index];
          // dmap_total_trilinear += tex3D<fp_type>(desolv_texture, coord_tex[0], coord_tex[1], coord_tex[2]) *
          //                         fabsf(l_ligand_charge[atom_index]);
          // const auto temp = tex3D<fp_type>(atom_textures[l_atom_tex_indexes[atom_index]],
          //                                        coord_tex[0],
          //                                        coord_tex[1],
          //                                        coord_tex[2]);
          //                                        emap_total_trilinear +=temp;
          elect_total_trilinear +=
              trilinear_interpolation_cuda(int_coord, electro_texture, coeffs) * ligand_charge[atom_index];
          dmap_total_trilinear += trilinear_interpolation_cuda(int_coord, desolv_texture, coeffs) *
                                  fabsf(ligand_charge[atom_index]);
          emap_total_trilinear +=
              trilinear_interpolation_cuda(int_coord, atom_textures[atom_tex_indexes[atom_index]], coeffs);
        }
      }
    }

    fp_type total_trilinear_eintcal = elect_total_trilinear + dmap_total_trilinear + emap_total_trilinear;
    __syncwarp();

    fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
    if (num_rotamers > 0)
      for (int nonbond_list = threadIdx.x; nonbond_list < ligand_num_nonbonds; nonbond_list += blockDim.x) {
        const int& a1 = ligand_nonbond_a1[nonbond_list];
        const int& a2 = ligand_nonbond_a2[nonbond_list];

        const auto diff_x                = ligand_x[a1] - ligand_x[a2];
        const auto diff_y                = ligand_y[a1] - ligand_y[a2];
        const auto diff_z                = ligand_z[a1] - ligand_z[a2];
        const fp_type distance_two       = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
        const fp_type distance_two_clamp = std::max(distance_two, RMIN_ELEC_SQUARE_CUDA);
        const fp_type distance           = sqrtf(distance_two_clamp);

        //  Calculate  Electrostatic  Energy
        const fp_type epsilon      = A + B / (fp_type{1} + rk * expf(lambda_B * distance));
        const fp_type r_dielectric = fp_type{1} / (distance * epsilon);
        const fp_type e_elec       = ligand_charge[a1] * ligand_charge[a2] * ELECSCALE *
                               autodock_parameters::coeff_estat * r_dielectric;
        elect_total_eintcal += e_elec;

        // Calcuate desolv
        const fp_type nb_desolv = (ligand_vol[a2] * (ligand_solpar[a1] + qsolpar * fabsf(ligand_charge[a1])) +
                                   ligand_vol[a1] * (ligand_solpar[a2] + qsolpar * fabsf(ligand_charge[a2])));

        const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                 expf(fp_type{-0.5} / (sigma_square_cuda) *distance_two_clamp) * nb_desolv;
        dmap_total_eintcal += e_desolv;

        fp_type e_vdW_Hb{0};
        if (distance_two_clamp < nbc2) {
          //   const int& hbond_i        = ligand_num_hbond[a1];
          //   const int& hbond_j        = ligand_num_hbond[a2];
          //   const fp_type& Rij_hb_i   = ligand_Rij_hb[a1];
          //   const fp_type& Rij_hb_j   = ligand_Rij_hb[a2];
          //   const fp_type& Rii_i      = ligand_Rii[a1];
          //   const fp_type& Rii_j      = ligand_Rii[a2];
          //   const fp_type& epsij_hb_i = ligand_epsij_hb[a1];
          //   const fp_type& epsij_hb_j = ligand_epsij_hb[a2];
          //   const fp_type& epsii_i    = ligand_epsii[a1];
          //   const fp_type& epsii_j    = ligand_epsii[a2];
          //
          //   int xA = 12;
          //   int xB = 6;
          //
          //   fp_type Rij{0}, epsij{0};
          //   if ((hbond_i == 1 || hbond_i == 2) && hbond_j > 2) {
          //     Rij   = Rij_hb_j;
          //     epsij = epsij_hb_j;
          //     xB    = 10;
          //   } else if ((hbond_i > 2) && (hbond_j == 1 || hbond_j == 2)) {
          //     Rij   = Rij_hb_i;
          //     epsij = epsij_hb_i;
          //     xB    = 10;
          //   } else {
          //     Rij   = (Rii_i + Rii_j) / fp_type{2};
          //     epsij = sqrtf(epsii_i * epsii_j);
          //   }
          //   if (xA != xB) {
          //     const fp_type tmp = epsij / (xA - xB);
          //     const fp_type cA  = tmp * powf(Rij, xA) * xB;
          //     const fp_type cB  = tmp * powf(Rij, xB) * xA;
          //
          //     const fp_type rA = powf(distance, static_cast<fp_type>(xA));
          //     const fp_type rB = powf(distance, static_cast<fp_type>(xB));
          //
          //     e_vdW_Hb = fminf(EINTCLAMP, (cA / rA - cB / rB));
          // }
          const int xA = xA_default;
          const int xB = ligand_nonbond_xB[nonbond_list];

          if (xA != xB) {
            const fp_type cA = ligand_nonbond_cA[nonbond_list];
            const fp_type cB = ligand_nonbond_cB[nonbond_list];

            const auto log_distance = std::log(distance);
            const fp_type rA        = std::exp(static_cast<fp_type>(xA) * log_distance);
            const fp_type rB        = std::exp(static_cast<fp_type>(xB) * log_distance);

            e_vdW_Hb = std::min(EINTCLAMP_CUDA, (cA / rA - cB / rB));
          }
        }
        emap_total_eintcal += e_vdW_Hb;
      }
    return elect_total_eintcal + emap_total_eintcal + dmap_total_eintcal + total_trilinear_eintcal;
  }
} // namespace mudock
