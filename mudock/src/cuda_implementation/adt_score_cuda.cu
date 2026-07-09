#include <array>
#include <cassert>
#include <memory>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/compute/adt_score_kernel.hpp>
#include <mudock/compute/devices_memory.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cuda_implementation/adt_score_cuda.cuh>
#include <mudock/cuda_implementation/cuda_texture.cuh>
#include <mudock/cuda_implementation/cuda_utils.cuh>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <stdexcept>
#include <stdio.h>

#define BUCKET_MULTIPLIER 3
#define FLATTENED_3D(x, y, z, index_x, index_xy) (index_xy * (z) + (y) * index_x + (x))

namespace mudock {
  __device__ static constexpr fp_type EINTCLAMP_CUDA{EINTCLAMP};
  __device__ static constexpr fp_type A{mehler_solmajer::A};
  __device__ static constexpr fp_type B{mehler_solmajer::B};
  __device__ static constexpr fp_type rk{mehler_solmajer::rk};
  __device__ static constexpr fp_type lambda_B{mehler_solmajer::lambda_B};

  __device__ static constexpr fp_type RMIN_ELEC_SQUARE_CUDA = RMIN_ELEC_SQUARE;
  __device__ static constexpr fp_type sigma_square_cuda     = sigma_square;

  __device__ __constant__ fp_type map_min_const[3];
  __device__ __constant__ fp_type map_max_const[3];
  __device__ __constant__ fp_type map_center_const[3];

  constexpr int k_max_devices = 16;

  device_memory_array<k_max_devices, cuda_texture_devices> cuda_texture_memory;
  device_memory_array<k_max_devices, fp_type> cuda_constant_memory;

  void init_device(const int dev,
                   const fp_type* map_min,
                   const fp_type* map_max,
                   const fp_type* map_center,
                   const int map_index_xyz,
                   const fp_type* map_grids) {
    // Thread-safe, exactly-once init per device:
    cuda_texture_memory.init(dev, map_index_xyz, num_autodock_grids(), map_grids);
    cuda_constant_memory.init(
        dev,
        std::function<std::unique_ptr<fp_type>()>([&]() {
          MUDOCK_CHECK(cudaSetDevice(dev));
          // MUDOCK_CHECK(cudaStreamCreate(&stream));
          MUDOCK_CHECK(
              cudaMemcpyToSymbol(map_min_const, map_min, 3 * sizeof(fp_type), 0, cudaMemcpyDeviceToDevice));
          MUDOCK_CHECK(
              cudaMemcpyToSymbol(map_max_const, map_max, 3 * sizeof(fp_type), 0, cudaMemcpyDeviceToDevice));
          MUDOCK_CHECK(cudaMemcpyToSymbol(map_center_const,
                                          map_center,
                                          3 * sizeof(fp_type),
                                          0,
                                          cudaMemcpyDeviceToDevice));
          MUDOCK_CHECK(cudaDeviceSynchronize());
          // TODO implement RAI to release constant memory
          return std::make_unique<fp_type>(0);
        }));
  }

  __device__ __forceinline__ fp_type dot3(const fp_type a[3], const fp_type b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }

  __device__ __forceinline__ void cross3(const fp_type a[3], const fp_type b[3], fp_type out[3]) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
  }

  __device__ __forceinline__ void get_grid_values(const fp_type *__restrict__ map,
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

  __device__ __forceinline__ fp_type trilinear_interpolation_cuda(const fp_type* __restrict__ map,
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
  __global__ void calc_energy(const int atom_stride,
                              const int scores_per_ligand,
                              const fp_type* __restrict__ scratch_x,
                              const fp_type* __restrict__ scratch_y,
                              const fp_type* __restrict__ scratch_z,
                              const fp_type* __restrict__ vol,
                              const fp_type* __restrict__ solpar,
                              const fp_type* __restrict__ charge,
                              const int* num_atoms_b,
                              const int* num_rotamers_b,
                              const int* num_nonbonds_b,
                              const int* __restrict__ nonbond_a1,
                              const int* __restrict__ nonbond_a2,
                              const fp_type* __restrict__ nonbond_cA,
                              const fp_type* __restrict__ nonbond_cB,
                              const int* __restrict__ nonbond_xB,
                              const int map_index_x,
                              const int map_index_xy,
                              const int map_index_xyz,
                              const fp_type* __restrict__ grid_maps,
                              const int* __restrict__ atom_tex_indexes,
                              fp_type* __restrict__ scores) {
    const fp_type* electro_map = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
    const fp_type* desolv_map  = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

    const int ligand_id       = blockIdx.x;
    const int local_thread_id = threadIdx.x;
    assert(blockDim.x == BLOCK_SIZE && warpSize == BLOCK_SIZE &&
           "Warpsize and the number of thread per block does not coincide");

    const int num_atoms    = num_atoms_b[ligand_id];
    const int num_nonbonds = num_nonbonds_b[ligand_id + 1] - num_nonbonds_b[ligand_id];
    const int num_rotamers = num_rotamers_b[ligand_id];
    const int stride       = ligand_id * atom_stride;

    const fp_type* l_scratch_x = scratch_x + stride * scores_per_ligand;
    const fp_type* l_scratch_y = scratch_y + stride * scores_per_ligand;
    const fp_type* l_scratch_z = scratch_z + stride * scores_per_ligand;
    const fp_type* l_vol       = vol + stride;
    const fp_type* l_solpar    = solpar + stride;
    const fp_type* l_charge    = charge + stride;
    // Point to the next population buffer
    const auto* l_atom_tex_indexes = atom_tex_indexes + stride;
    const int* l_nonbond_a1        = nonbond_a1 + num_nonbonds_b[ligand_id];
    const int* l_nonbond_a2        = nonbond_a2 + num_nonbonds_b[ligand_id];
    const fp_type* l_nonbond_cA    = nonbond_cA + num_nonbonds_b[ligand_id];
    const fp_type* l_nonbond_cB    = nonbond_cB + num_nonbonds_b[ligand_id];
    const int* l_nonbond_xB        = nonbond_xB + num_nonbonds_b[ligand_id];

    fp_type* scores_l = scores + ligand_id * scores_per_ligand;

    for (int scores_index = 0; scores_index < scores_per_ligand; ++scores_index) {
      // Copy original coordinates
      const fp_type* ligand_x = l_scratch_x + scores_index * atom_stride;
      const fp_type* ligand_y = l_scratch_y + scores_index * atom_stride;
      const fp_type* ligand_z = l_scratch_z + scores_index * atom_stride;

      // Calculate energy
      fp_type elect_total_trilinear = 0, emap_total_trilinear = 0, dmap_total_trilinear = 0;
      MUDOCK_PRAGMA_UNROLL(MUDOCK_ATOM_LOOP_UNROLL_FACTOR(MAX_ATOMS, BLOCK_SIZE))
      for (int atom_index = threadIdx.x; atom_index < MAX_ATOMS; atom_index += BLOCK_SIZE) {
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
            coord_tex[0]       = (coord_tex[0] - map_min_const[0]) * inv_spacing,
            coord_tex[1]       = (coord_tex[1] - map_min_const[1]) * inv_spacing;
            coord_tex[2]       = (coord_tex[2] - map_min_const[2]) * inv_spacing;
            const auto& charge = l_charge[atom_index];
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
            const int base_index    = FLATTENED_3D(u0, v0, w0, map_index_x, map_index_xy);
            const fp_type* atom_map = grid_maps + l_atom_tex_indexes[atom_index];

            elect_total_trilinear +=
                trilinear_interpolation_cuda(electro_map + base_index, coeffs, map_index_x, map_index_xy) *
                charge;
            dmap_total_trilinear +=
                trilinear_interpolation_cuda(desolv_map + base_index, coeffs, map_index_x, map_index_xy) *
                fabsf(charge);
            emap_total_trilinear +=
                trilinear_interpolation_cuda(atom_map + base_index, coeffs, map_index_x, map_index_xy);
          }
        }
      }

      fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
      if (num_rotamers > 0)
        for (int nonbond_list = threadIdx.x; nonbond_list < num_nonbonds; nonbond_list += blockDim.x) {
          const int& a1    = l_nonbond_a1[nonbond_list];
          const int& a2    = l_nonbond_a2[nonbond_list];
          const auto& a1_c = l_charge[a1];
          const auto& a2_c = l_charge[a2];

          const auto diff_x          = ligand_x[a1] - ligand_x[a2];
          const auto diff_y          = ligand_y[a1] - ligand_y[a2];
          const auto diff_z          = ligand_z[a1] - ligand_z[a2];
          const fp_type distance_two = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
          const fp_type distance_two_clamp =
              distance_two > RMIN_ELEC_SQUARE_CUDA ? distance_two : RMIN_ELEC_SQUARE_CUDA;
          const fp_type distance = sqrtf(distance_two_clamp);

          //  Calculate  Electrostatic  Energy
          const fp_type epsilon      = A + B / (fp_type{1} + rk * expf(lambda_B * distance));
          const fp_type r_dielectric = fp_type{1} / (distance * epsilon);
          const fp_type e_elec = a1_c * a2_c * ELECSCALE * autodock_parameters::coeff_estat * r_dielectric;
          elect_total_eintcal += e_elec;

          // Calcuate desolv
          const fp_type nb_desolv = (l_vol[a2] * (l_solpar[a1] + qsolpar * fabsf(a1_c)) +
                                     l_vol[a1] * (l_solpar[a2] + qsolpar * fabsf(a2_c)));

          const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                   expf(fp_type{-0.5} / (sigma_square_cuda) *distance_two_clamp) * nb_desolv;
          dmap_total_eintcal += e_desolv;

          fp_type e_vdW_Hb{0};
          if (distance_two_clamp < nbc2) {
            const int xA = xA_default;
            const int xB = l_nonbond_xB[nonbond_list];

            if (xA != xB) {
              const fp_type cA = l_nonbond_cA[nonbond_list];
              const fp_type cB = l_nonbond_cB[nonbond_list];

              const auto log_distance = logf(distance);
              const fp_type rA        = expf(static_cast<fp_type>(xA) * log_distance);
              const fp_type rB        = expf(static_cast<fp_type>(xB) * log_distance);

              e_vdW_Hb = (cA / rA - cB / rB);
              e_vdW_Hb = EINTCLAMP_CUDA < e_vdW_Hb ? EINTCLAMP_CUDA : e_vdW_Hb;
            }
          }
          emap_total_eintcal += e_vdW_Hb;
        }

      fp_type total_energy = emap_total_eintcal + elect_total_eintcal + dmap_total_eintcal +
                             emap_total_trilinear + elect_total_trilinear + dmap_total_trilinear;

      MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
      for (int offset = BLOCK_SIZE / 2; offset > 0; offset /= 2) {
        total_energy += __shfl_down_sync(0xffffffff, total_energy, offset);
      }

      if (local_thread_id == 0) {
        const fp_type tors_free_energy = num_rotamers * autodock_parameters::coeff_tors;
        scores_l[scores_index]         = total_energy + tors_free_energy;
      }
    }
  }


template<int MAX_ATOMS>
__global__ void calc_gradient(const int batch_atoms,
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
                              const int *__restrict__ map_offsets_b,
                              const int map_index_x,
                              const int map_index_xy,
                              const int map_index_xyz,
                              gradient *__restrict__ gradients_b,
                              int *__restrict__ active_individuals_b) {
    const fp_type *electro_map = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
    const fp_type *desolv_map  = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

    const int ligand_id       = blockIdx.x;
    assert(blockDim.x == BLOCK_SIZE && warpSize == BLOCK_SIZE &&
           "Warpsize and the number of thread per block does not coincide");

    const int atom_stride                  = ligand_id * batch_atoms;
    const int num_atoms                    = num_atoms_b[ligand_id];
    const int num_nonbonds                 = num_nonbonds_b[ligand_id + 1] - num_nonbonds_b[ligand_id];
    const int num_rotamers                 = num_rotamers_b[ligand_id];
    const int batch_rotamers               = batch_atoms - 3; // TODO L why is this the max number of rotamers?
    const fp_type *__restrict__ scratch_x  = x_scratch_b + atom_stride * individuals_per_ligand;
    const fp_type *__restrict__ scratch_y  = y_scratch_b + atom_stride * individuals_per_ligand;
    const fp_type *__restrict__ scratch_z  = z_scratch_b + atom_stride * individuals_per_ligand;
    const chromosome* __restrict__ ligand_chromosomes = chromosomes_b + ligand_id * individuals_per_ligand;
    const fp_type *__restrict__ vol_l      = vols_b + atom_stride;
    const fp_type *__restrict__ solpar_l   = solpars_b + atom_stride;
    const fp_type *__restrict__ charge_l   = charges_b + atom_stride;
    const int *__restrict__ map_offsets_l  = map_offsets_b + atom_stride;
    const int *__restrict__ nonbond_a1_l   = nonbond_a1_b + num_nonbonds_b[ligand_id];
    const int *__restrict__ nonbond_a2_l   = nonbond_a2_b + num_nonbonds_b[ligand_id];
    const fp_type *nonbond_cA_l            = nonbond_cA_b + num_nonbonds_b[ligand_id];
    const fp_type *nonbond_cB_l            = nonbond_cB_b + num_nonbonds_b[ligand_id];
    const int *nonbond_xB_l                = nonbond_xB_b + num_nonbonds_b[ligand_id];
    const int* fragments                   = ligand_fragments_b + ligand_fragments_start_b[ligand_id];
    const int* frag_start_indices          = frag_start_indices_b + frag_indices_start_b[ligand_id];
    const int* frag_stop_indices           = frag_stop_indices_b + frag_indices_start_b[ligand_id];
    gradient *__restrict__ gradients_l     = gradients_b + ligand_id * individuals_per_ligand;
    int *__restrict__ active_individuals_l = active_individuals_b + ligand_id * individuals_per_ligand;
    

    // TODO L very Important: fix this magic number with actual values
    fp_type dE_dX[3 * MAX_ATOMS];
    fp_type grad[6 + MAX_ATOMS];    // gradient = dE/dx, dE/dy, dE/dz, dE/dalpha, dE/dbeta, dE/dgamma, dE/d_tors_1, ..., dE/d_tors_n

    // TODO L now each block compute a ligand and each thread in the block compute a set of individuals, strided by blockDim
    // and no parallelization of atoms. Make some tests and see if it is better to divide like in calc_energy, where all threads
    // work on all population, but divide computation on atoms
    for (int individual_index = threadIdx.x; individual_index < individuals_per_ligand; individual_index += blockDim.x) {
      
      if (!active_individuals_l[individual_index]){
        continue;
      }

      // Reinitialize to 0s all gradient components vectors
      for (int i = 0; i < 3 * batch_atoms; ++i) dE_dX[i] = 0;
      for (int i = 0; i < 6 + batch_rotamers; ++i) grad[i] = 0;

      const fp_type *__restrict__ scratch_x_l = scratch_x + individual_index * batch_atoms;
      const fp_type *__restrict__ scratch_y_l = scratch_y + individual_index * batch_atoms;
      const fp_type *__restrict__ scratch_z_l = scratch_z + individual_index * batch_atoms;

      for (int index = 0; index < num_atoms; ++index) {
        fp_type coord[3]{scratch_x_l[index], scratch_y_l[index], scratch_z_l[index]};

        if (coord[0] < map_min_const[0] || coord[0] > map_max_const[0] || coord[1] < map_min_const[1] ||
          coord[1] > map_max_const[1] || coord[2] < map_min_const[2] || coord[2] > map_max_const[2]) {
          const auto diff_x = coord[0] - map_center_const[0];
          const auto diff_y = coord[1] - map_center_const[1];
          const auto diff_z = coord[2] - map_center_const[2];
          const fp_type penalty_factor = 2 * 2 * ENERGYPENALTY;

          dE_dX[3*index]     += penalty_factor * diff_x;
          dE_dX[3*index + 1] += penalty_factor * diff_y;
          dE_dX[3*index + 2] += penalty_factor * diff_z;
        } else {
          const auto &atom_charge = charge_l[index];
          const fp_type *atom_map = grid_maps + map_offsets_l[index];

          coord[0] = (coord[0] - map_min_const[0]) * inv_spacing;
          coord[1] = (coord[1] - map_min_const[1]) * inv_spacing;
          coord[2] = (coord[2] - map_min_const[2]) * inv_spacing;

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
          dE_dX[3*index]     += constant_factor_el * (p1w * (p1v * (grid_values[4] - grid_values[0]) + p0v * (grid_values[6] - grid_values[2])) + p0w * (p1v * (grid_values[5] - grid_values[1]) + p0v * (grid_values[7] - grid_values[3])));
          dE_dX[3*index + 1] += constant_factor_el * (p1w * (p1u * (grid_values[2] - grid_values[0]) + p0u * (grid_values[6] - grid_values[4])) + p0w * (p1u * (grid_values[3] - grid_values[1]) + p0u * (grid_values[7] - grid_values[5])));
          dE_dX[3*index + 2] += constant_factor_el * (p1v * (p1u * (grid_values[1] - grid_values[0]) + p0u * (grid_values[5] - grid_values[4])) + p0v * (p1u * (grid_values[3] - grid_values[2]) + p0u * (grid_values[7] - grid_values[6])));

          get_grid_values(atom_map + base_index, map_index_x, map_index_xy, grid_values);
          dE_dX[3*index]     += inv_spacing * (p1w * (p1v * (grid_values[4] - grid_values[0]) + p0v * (grid_values[6] - grid_values[2])) + p0w * (p1v * (grid_values[5] - grid_values[1]) + p0v * (grid_values[7] - grid_values[3])));
          dE_dX[3*index + 1] += inv_spacing * (p1w * (p1u * (grid_values[2] - grid_values[0]) + p0u * (grid_values[6] - grid_values[4])) + p0w * (p1u * (grid_values[3] - grid_values[1]) + p0u * (grid_values[7] - grid_values[5])));
          dE_dX[3*index + 2] += inv_spacing * (p1v * (p1u * (grid_values[1] - grid_values[0]) + p0u * (grid_values[5] - grid_values[4])) + p0v * (p1u * (grid_values[3] - grid_values[2]) + p0u * (grid_values[7] - grid_values[6])));            
          
          get_grid_values(desolv_map + base_index, map_index_x, map_index_xy, grid_values);
          const fp_type constant_factor_des = inv_spacing * std::fabs(atom_charge);
          dE_dX[3*index]     += constant_factor_des * (p1w * (p1v * (grid_values[4] - grid_values[0]) + p0v * (grid_values[6] - grid_values[2])) + p0w * (p1v * (grid_values[5] - grid_values[1]) + p0v * (grid_values[7] - grid_values[3])));
          dE_dX[3*index + 1] += constant_factor_des * (p1w * (p1u * (grid_values[2] - grid_values[0]) + p0u * (grid_values[6] - grid_values[4])) + p0w * (p1u * (grid_values[3] - grid_values[1]) + p0u * (grid_values[7] - grid_values[5])));
          dE_dX[3*index + 2] += constant_factor_des * (p1v * (p1u * (grid_values[1] - grid_values[0]) + p0u * (grid_values[5] - grid_values[4])) + p0v * (p1u * (grid_values[3] - grid_values[2]) + p0u * (grid_values[7] - grid_values[6])));
        }
      }
      
      if (num_rotamers > 0) {
        for (int i = 0; i < num_nonbonds; ++i) {
          const int &a1 = nonbond_a1_l[i];
          const int &a2 = nonbond_a2_l[i];

          const auto diff_x                = scratch_x_l[a1] - scratch_x_l[a2];
          const auto diff_y                = scratch_y_l[a1] - scratch_y_l[a2];
          const auto diff_z                = scratch_z_l[a1] - scratch_z_l[a2];
          const fp_type distance_two       = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
          const fp_type distance_two_clamp = std::max(distance_two, RMIN_ELEC_SQUARE_CUDA);
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

          dE_dX[3*a1] += gx;
          dE_dX[3*a1 + 1] += gy;
          dE_dX[3*a1 + 2] += gz;

          dE_dX[3*a2] -= gx;
          dE_dX[3*a2 + 1] -= gy;
          dE_dX[3*a2 + 2] -= gz;

        }
      }


      fp_type ligand_COM[3] = {0};
      for (int i = 0; i < num_atoms; ++i) {
        ligand_COM[0] += scratch_x_l[i];
        ligand_COM[1] += scratch_y_l[i];
        ligand_COM[2] += scratch_z_l[i];
      }
      ligand_COM[0] /= static_cast<fp_type>(num_atoms);
      ligand_COM[1] /= static_cast<fp_type>(num_atoms);
      ligand_COM[2] /= static_cast<fp_type>(num_atoms);

      fp_type tau[3] = {0};
      for (int i = 0; i < num_atoms; ++i) {
        fp_type r[3] = {0};
        r[0] = scratch_x_l[i] - ligand_COM[0];
        r[1] = scratch_y_l[i] - ligand_COM[1];
        r[2] = scratch_z_l[i] - ligand_COM[2];
        fp_type torque[3] = {0};
        cross3(r, &dE_dX[3*i], torque);
        tau[0] += torque[0];
        tau[1] += torque[1];
        tau[2] += torque[2];
      }

      const chromosome &chrom = ligand_chromosomes[individual_index];
      const fp_type rad_beta = deg_to_rad(chrom[4]);
      const fp_type rad_gamma = deg_to_rad(chrom[5]);

      const fp_type sin_beta = std::sin(rad_beta);
      const fp_type cos_beta = std::cos(rad_beta);
      const fp_type sin_gamma = std::sin(rad_gamma);
      const fp_type cos_gamma = std::cos(rad_gamma);

      fp_type u_gamma[3];
      u_gamma[0] = fp_type{0};
      u_gamma[1] = fp_type{0};
      u_gamma[2] = fp_type{1};

      fp_type u_beta[3];
      u_beta[0] = -sin_gamma;
      u_beta[1] = cos_gamma;
      u_beta[2] = fp_type{0};

      fp_type u_alpha[3];
      u_alpha[0] = cos_beta * cos_gamma;
      u_alpha[1] = cos_beta * sin_gamma;
      u_alpha[2] = -sin_beta;

      grad[3] = dot3(tau, u_alpha);
      grad[4] = dot3(tau, u_beta);
      grad[5] = dot3(tau, u_gamma);

      for (int index = 0; index < num_atoms; ++index) {
        grad[0] += dE_dX[3*index];
        grad[1] += dE_dX[3*index + 1];
        grad[2] += dE_dX[3*index + 2];
      }

      // Torsion derivatives
      for (int t = 0; t < num_rotamers; ++t) {
        const int* frag_mask = fragments + t * num_atoms;

        const int a1 = frag_start_indices[t];
        const int a2 = frag_stop_indices[t];

        fp_type axis[3];
        axis[0] = scratch_x_l[a2] - scratch_x_l[a1];
        axis[1] = scratch_y_l[a2] - scratch_y_l[a1];
        axis[2] = scratch_z_l[a2] - scratch_z_l[a1];

        const fp_type norm = std::sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);

        const fp_type inv_norm = fp_type{1} / norm;
        axis[0] *= inv_norm;
        axis[1] *= inv_norm;
        axis[2] *= inv_norm;

        fp_type r[3] = {0};
        fp_type temp[3] = {0};
        for (int i = 0; i < num_atoms; ++i) {
          if (frag_mask[i] != 0) {
            r[0] = scratch_x_l[i] - scratch_x_l[a1];
            r[1] = scratch_y_l[i] - scratch_y_l[a1];
            r[2] = scratch_z_l[i] - scratch_z_l[a1];
            cross3(axis, r, temp);
            grad[6 + t] += dot3(&dE_dX[3*i], temp);
          }
        }
      }

      // Write gradient on the buffer
      gradient &grad_l = gradients_l[individual_index];
      for (int i = 0; i < 6 + num_rotamers; ++i) {
        grad_l[i] = grad[i];
      } 
    }
  };

  template<>
  void adt_score_kernel<queue_cuda>::operator()() {
    const int dev_id = q->get_id();
    init_device(dev_id, minimum, maximum, center, map_index_xyz, grid_maps);

    void* args[] = {(void*) &batch_atoms,
                    (void*) &scores_per_ligand,
                    (void*) &x_scratch_b,
                    (void*) &y_scratch_b,
                    (void*) &z_scratch_b,
                    (void*) &vols_b,
                    (void*) &solpars_b,
                    (void*) &charges_b,
                    (void*) &num_atoms_b,
                    (void*) &num_rotamers_b,
                    (void*) &num_nonbonds_b,
                    (void*) &nonbond_a1_b,
                    (void*) &nonbond_a2_b,
                    (void*) &nonbond_cA_b,
                    (void*) &nonbond_cB_b,
                    (void*) &nonbond_xB_b,
                    (void*) &map_index_x,
                    (void*) &map_index_xy,
                    (void*) &map_index_xyz,
                    (void*) &(*cuda_texture_memory.v[dev_id].data).tex_dev,
                    (void*) &map_offsets_b,
                    (void*) &scores_b};
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->launch_kernel((void*) calc_energy<max_atoms>, args, batch_ligands, BLOCK_SIZE);
        },
        batch_atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());

  }; // namespace mudock

  template<>
  void adt_gradient_kernel<queue_cuda>::operator()() {
    const int dev_id = q->get_id();
    // init_device(dev_id, minimum, maximum, center, map_index_xyz, grid_maps);
    
    void* args[] = {(void*) &batch_atoms,
                    (void*) &batch_ligands,
                    (void*) &scores_per_ligand,
                    (void*) &x_scratch_b,
                    (void*) &y_scratch_b,
                    (void*) &z_scratch_b,
                    (void*) &chromosomes_b,
                    (void*) &ligand_fragments_b,
                    (void*) &ligand_fragments_start_b,
                    (void*) &frag_indices_start_b,
                    (void*) &frag_start_indices_b,
                    (void*) &frag_stop_indices_b,
                    (void*) &vols_b,
                    (void*) &solpars_b,
                    (void*) &charges_b,
                    (void*) &num_atoms_b,
                    (void*) &num_rotamers_b,
                    (void*) &num_nonbonds_b,
                    (void*) &nonbond_a1_b,
                    (void*) &nonbond_a2_b,
                    (void*) &nonbond_cA_b,
                    (void*) &nonbond_cB_b,
                    (void*) &nonbond_xB_b,
                    (void*) &(*cuda_texture_memory.v[dev_id].data).tex_dev,
                    (void*) &map_offsets_b,
                    (void*) &map_index_x,
                    (void*) &map_index_xy,
                    (void*) &map_index_xyz,
                    (void*) &gradients_b,
                    (void*) &active_individuals_b};
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->launch_kernel((void*) calc_gradient<max_atoms>, args, batch_ligands, BLOCK_SIZE);
        },
        batch_atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

  template<int MAX_ATOMS>
  batch_multiple get_evaluate_fitness_batch(const int device_id) {
    return get_kernel_batch_multiple_cuda<calc_energy<MAX_ATOMS>>(device_id,
                                                                  BLOCK_SIZE,
                                                                  0,
                                                                  "adt_score::calc_energy");
  }

  template<>
  batch_multiple get_adt_score_batch_multiple<queue_cuda>(const int atoms, std::shared_ptr<queue_cuda> q_b) {
    batch_multiple bucket_multiple{};
    const int device_id = q_b->get_id();
    constexpr_for<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>([&](const auto atom_index) {
      const auto n_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
      if (atoms == n_atoms)
        bucket_multiple = get_evaluate_fitness_batch<n_atoms>(device_id);
    });
    if (bucket_multiple.total_multiple() <= 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");
    return normalize_batch_multiple(bucket_multiple);
  }

} // namespace mudock
