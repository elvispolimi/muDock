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

  __device__ __forceinline__ fp_type warp_allreduce_sum(fp_type value) {
    MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
    for (int offset = BLOCK_SIZE / 2; offset > 0; offset /= 2) {
      value += __shfl_xor_sync(0xffffffff, value, offset);
    }
    return value;
  }

  // One full warp owns the shared gradient arrays. All 32 lanes must call
  // this helper, using atom == -1 for an inactive pair. Within each phase,
  // only one leader writes each atom; the final barrier orders subsequent
  // phases, so plain shared-memory additions are sufficient.
  // Requires compute capability >= 7.0 (__match_any_sync).
  __device__ __forceinline__ void accumulate_without_atomics(
      fp_type* dx, fp_type* dy, fp_type* dz,
      const int atom, const fp_type gx, const fp_type gy, const fp_type gz) {
    constexpr unsigned full_mask = 0xffffffffu;
    const int lane = threadIdx.x;
    const unsigned peers = __match_any_sync(full_mask, atom);

    if (atom >= 0) {
      const int leader = __ffs(peers) - 1;
      fp_type sx = gx;
      fp_type sy = gy;
      fp_type sz = gz;

      if (__popc(peers) > 1) {
        sx = sy = sz = fp_type{0};
        unsigned remaining = peers;
        while (remaining != 0) {
          const int source = __ffs(remaining) - 1;
          sx += __shfl_sync(peers, gx, source);
          sy += __shfl_sync(peers, gy, source);
          sz += __shfl_sync(peers, gz, source);
          remaining &= remaining - 1;
        }
      }

      if (lane == leader) {
        dx[atom] += sx;
        dy[atom] += sy;
        dz[atom] += sz;
      }
    }

    // Every lane participates, including lanes with no valid pair. This
    // separates endpoint updates and prevents overlap between loop rounds.
    __syncwarp(full_mask);
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
                              fp_type* __restrict__ scores,
                              int *__restrict__ converged_ligands_b) {
    const fp_type* electro_map = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
    const fp_type* desolv_map  = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

    const int ligand_id       = blockIdx.x;
    const int local_thread_id = threadIdx.x;
    assert(blockDim.x == BLOCK_SIZE && warpSize == BLOCK_SIZE &&
           "Warpsize and the number of thread per block does not coincide");

    if (converged_ligands_b[ligand_id]) return;

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

          if (coord_tex[0] < map_min_const[0] || coord_tex[0] >= map_max_const[0] ||
              coord_tex[1] < map_min_const[1] || coord_tex[1] >= map_max_const[1] ||
              coord_tex[2] < map_min_const[2] || coord_tex[2] >= map_max_const[2]) {
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
                                int *__restrict__ ls_selected_b,
                                int *__restrict__ converged_ligands_b) {
    const fp_type *electro_map = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
    const fp_type *desolv_map  = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

    const int ligand_id       = blockIdx.x;
    const int local_thread_id = threadIdx.x;
    assert(blockDim.x == BLOCK_SIZE && warpSize == BLOCK_SIZE &&
           "Warpsize and the number of thread per block does not coincide");

    // Block-uniform early exit: all lanes leave together
    if (converged_ligands_b[ligand_id]) return;

    const int atom_stride                  = ligand_id * batch_atoms;
    const int num_atoms                    = num_atoms_b[ligand_id];
    const int num_nonbonds                 = num_nonbonds_b[ligand_id + 1] - num_nonbonds_b[ligand_id];
    const int num_rotamers                 = num_rotamers_b[ligand_id];
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
    int *__restrict__ ls_selected_l        = ls_selected_b + ligand_id * individuals_per_ligand;

    const fp_type angle_scale = deg_to_rad(fp_type{1});
    const fp_type inv_n       = fp_type{1} / static_cast<fp_type>(num_atoms);

    // Per-atom gradient, shared by the 32 lanes of the block
    __shared__ fp_type dE_dx[MAX_ATOMS];
    __shared__ fp_type dE_dy[MAX_ATOMS];
    __shared__ fp_type dE_dz[MAX_ATOMS];

    // All lanes work together on one individual at a time
    for (int individual_index = 0; individual_index < individuals_per_ligand; ++individual_index) {
      if (!ls_selected_l[individual_index]) continue; // block-uniform

      const fp_type *__restrict__ scratch_x_l = scratch_x + individual_index * batch_atoms;
      const fp_type *__restrict__ scratch_y_l = scratch_y + individual_index * batch_atoms;
      const fp_type *__restrict__ scratch_z_l = scratch_z + individual_index * batch_atoms;

      // ---------------------------------------------------------------------------------------------
      // 1. Grid (intermolecular) contribution: one thread per atom, writes only its own dE entries
      // ---------------------------------------------------------------------------------------------
      for (int index = local_thread_id; index < num_atoms; index += BLOCK_SIZE) {
        fp_type gx{0}, gy{0}, gz{0};
        fp_type coord[3]{scratch_x_l[index], scratch_y_l[index], scratch_z_l[index]};

        if (coord[0] < map_min_const[0] || coord[0] >= map_max_const[0] ||
            coord[1] < map_min_const[1] || coord[1] >= map_max_const[1] ||
            coord[2] < map_min_const[2] || coord[2] >= map_max_const[2]) {
          const auto diff_x = coord[0] - map_center_const[0];
          const auto diff_y = coord[1] - map_center_const[1];
          const auto diff_z = coord[2] - map_center_const[2];
          const fp_type penalty_factor = 2 * 2 * ENERGYPENALTY;

          gx += penalty_factor * diff_x;
          gy += penalty_factor * diff_y;
          gz += penalty_factor * diff_z;
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

          const int base_index = FLATTENED_3D(u0, v0, w0, map_index_x, map_index_xy);
          fp_type grid_values[8];

          get_grid_values(electro_map + base_index, map_index_x, map_index_xy, grid_values);
          const fp_type constant_factor_el = inv_spacing * atom_charge;
          gx += constant_factor_el * (p1w * (p1v * (grid_values[4] - grid_values[0]) + p0v * (grid_values[6] - grid_values[2])) + p0w * (p1v * (grid_values[5] - grid_values[1]) + p0v * (grid_values[7] - grid_values[3])));
          gy += constant_factor_el * (p1w * (p1u * (grid_values[2] - grid_values[0]) + p0u * (grid_values[6] - grid_values[4])) + p0w * (p1u * (grid_values[3] - grid_values[1]) + p0u * (grid_values[7] - grid_values[5])));
          gz += constant_factor_el * (p1v * (p1u * (grid_values[1] - grid_values[0]) + p0u * (grid_values[5] - grid_values[4])) + p0v * (p1u * (grid_values[3] - grid_values[2]) + p0u * (grid_values[7] - grid_values[6])));

          get_grid_values(atom_map + base_index, map_index_x, map_index_xy, grid_values);
          gx += inv_spacing * (p1w * (p1v * (grid_values[4] - grid_values[0]) + p0v * (grid_values[6] - grid_values[2])) + p0w * (p1v * (grid_values[5] - grid_values[1]) + p0v * (grid_values[7] - grid_values[3])));
          gy += inv_spacing * (p1w * (p1u * (grid_values[2] - grid_values[0]) + p0u * (grid_values[6] - grid_values[4])) + p0w * (p1u * (grid_values[3] - grid_values[1]) + p0u * (grid_values[7] - grid_values[5])));
          gz += inv_spacing * (p1v * (p1u * (grid_values[1] - grid_values[0]) + p0u * (grid_values[5] - grid_values[4])) + p0v * (p1u * (grid_values[3] - grid_values[2]) + p0u * (grid_values[7] - grid_values[6])));

          get_grid_values(desolv_map + base_index, map_index_x, map_index_xy, grid_values);
          const fp_type constant_factor_des = inv_spacing * std::fabs(atom_charge);
          gx += constant_factor_des * (p1w * (p1v * (grid_values[4] - grid_values[0]) + p0v * (grid_values[6] - grid_values[2])) + p0w * (p1v * (grid_values[5] - grid_values[1]) + p0v * (grid_values[7] - grid_values[3])));
          gy += constant_factor_des * (p1w * (p1u * (grid_values[2] - grid_values[0]) + p0u * (grid_values[6] - grid_values[4])) + p0w * (p1u * (grid_values[3] - grid_values[1]) + p0u * (grid_values[7] - grid_values[5])));
          gz += constant_factor_des * (p1v * (p1u * (grid_values[1] - grid_values[0]) + p0u * (grid_values[5] - grid_values[4])) + p0v * (p1u * (grid_values[3] - grid_values[2]) + p0u * (grid_values[7] - grid_values[6])));
        }
        // plain store also initializes the entry (replaces the explicit zeroing), no race: own atom
        dE_dx[index] = gx;
        dE_dy[index] = gy;
        dE_dz[index] = gz;
      }
      // Make the grid contributions visible before nonbonded accumulation.
      __syncwarp();

      // ---------------------------------------------------------------------------------------------
      // 2. Intramolecular non-bonded contribution: same pair-to-lane mapping,
      //    with atom-grouped updates and no atomic operations. The loop bounds
      //    are warp-uniform; inactive lanes still execute both helper calls.
      // ---------------------------------------------------------------------------------------------
      if (num_rotamers > 0) {
        for (int base = 0; base < num_nonbonds; base += BLOCK_SIZE) {
          const int i = base + local_thread_id;
          int target_a1 = -1;
          int target_a2 = -1;
          fp_type gx{0}, gy{0}, gz{0};

          if (i < num_nonbonds) {
            const int a1 = nonbond_a1_l[i];
            const int a2 = nonbond_a2_l[i];

            const auto diff_x = scratch_x_l[a1] - scratch_x_l[a2];
            const auto diff_y = scratch_y_l[a1] - scratch_y_l[a2];
            const auto diff_z = scratch_z_l[a1] - scratch_z_l[a2];
            const fp_type distance_two = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;

            // Preserve the original skip condition, but never skip a warp
            // collective. Invalid pairs retain the -1 endpoint sentinels.
            if (!(distance_two < RMIN_ELEC_SQUARE)) {
              target_a1 = a1;
              target_a2 = a2;

              const fp_type distance_two_clamp = distance_two;
              const fp_type distance           = std::sqrt(distance_two_clamp);
              const fp_type inv_r = fp_type{1} / distance;

              const fp_type dir_x = diff_x * inv_r;
              const fp_type dir_y = diff_y * inv_r;
              const fp_type dir_z = diff_z * inv_r;

              const fp_type exp_term = std::exp(mehler_solmajer::lambda_B * distance);
              const fp_type denom = (fp_type{1} + mehler_solmajer::rk * exp_term);
              const fp_type epsilon = mehler_solmajer::A + mehler_solmajer::B / denom;
              const fp_type epsilon_prime = -mehler_solmajer::B * mehler_solmajer::rk * mehler_solmajer::lambda_B * exp_term / (denom * denom);
              const fp_type C = charge_l[a1] * charge_l[a2] * ELECSCALE * autodock_parameters::coeff_estat;

              const fp_type f = distance * epsilon;
              const fp_type f_prime = epsilon + distance * epsilon_prime;
              const fp_type dE_dr_elec = -C * f_prime / (f * f);

              const fp_type nb_desolv = (vol_l[a2] * (solpar_l[a1] + qsolpar * std::fabs(charge_l[a1])) +
                                         vol_l[a1] * (solpar_l[a2] + qsolpar * std::fabs(charge_l[a2])));
              const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                       std::exp(fp_type{-0.5} / (sigma_square) * distance_two_clamp) * nb_desolv;
              const fp_type dE_dr_desolv = (-distance / sigma_square) * e_desolv;

              fp_type dE_dr_vdw = 0;
              if (distance_two_clamp < nbc2) {
                const int xA = xA_default;
                const int xB = nonbond_xB_l[i];

                if (xA != xB) {
                  const fp_type cA = nonbond_cA_l[i];
                  const fp_type cB = nonbond_cB_l[i];

                  const auto log_distance = std::log(distance);
                  const fp_type rA        = std::exp(-static_cast<fp_type>(xA) * log_distance);
                  const fp_type rB        = std::exp(-static_cast<fp_type>(xB) * log_distance);

                  const fp_type rA1 = rA * inv_r;
                  const fp_type rB1 = rB * inv_r;
                  const fp_type e_vdW_Hb_unclamped = cA * rA - cB * rB;
                  if (e_vdW_Hb_unclamped < EINTCLAMP) {
                    dE_dr_vdw = -static_cast<fp_type>(xA) * cA * rA1 + static_cast<fp_type>(xB) * cB * rB1;
                  }
                }
              }

              const fp_type dE_dr = dE_dr_elec + dE_dr_desolv + dE_dr_vdw;
              gx = dE_dr * dir_x;
              gy = dE_dr * dir_y;
              gz = dE_dr * dir_z;
            }
          }

          // Each phase has one writer per target atom. The helper ends with
          // a full-warp barrier, ordering a1, a2 and the next pair chunk.
          accumulate_without_atomics(dE_dx, dE_dy, dE_dz, target_a1, gx, gy, gz);
          accumulate_without_atomics(dE_dx, dE_dy, dE_dz, target_a2, -gx, -gy, -gz);
        }
      }
      // Ensure every final per-atom gradient is visible to the next phase.
      __syncwarp();

      // ---------------------------------------------------------------------------------------------
      // 3. Center of mass, torque, translation gradient: partial sums per lane + warp all-reduce
      // ---------------------------------------------------------------------------------------------
      fp_type com_x{0}, com_y{0}, com_z{0};
      fp_type gt_x{0}, gt_y{0}, gt_z{0};
      for (int i = local_thread_id; i < num_atoms; i += BLOCK_SIZE) {
        com_x += scratch_x_l[i];
        com_y += scratch_y_l[i];
        com_z += scratch_z_l[i];
        gt_x += dE_dx[i];
        gt_y += dE_dy[i];
        gt_z += dE_dz[i];
      }
      com_x = warp_allreduce_sum(com_x) * inv_n;
      com_y = warp_allreduce_sum(com_y) * inv_n;
      com_z = warp_allreduce_sum(com_z) * inv_n;
      const fp_type g_trasl_x = warp_allreduce_sum(gt_x);
      const fp_type g_trasl_y = warp_allreduce_sum(gt_y);
      const fp_type g_trasl_z = warp_allreduce_sum(gt_z);

      fp_type tau_x{0}, tau_y{0}, tau_z{0};
      for (int i = local_thread_id; i < num_atoms; i += BLOCK_SIZE) {
        const fp_type rx = scratch_x_l[i] - com_x;
        const fp_type ry = scratch_y_l[i] - com_y;
        const fp_type rz = scratch_z_l[i] - com_z;
        // Cross product: torque = r x dE_dX
        tau_x += ry * dE_dz[i] - rz * dE_dy[i];
        tau_y += rz * dE_dx[i] - rx * dE_dz[i];
        tau_z += rx * dE_dy[i] - ry * dE_dx[i];
      }
      fp_type tau[3];
      tau[0] = warp_allreduce_sum(tau_x);
      tau[1] = warp_allreduce_sum(tau_y);
      tau[2] = warp_allreduce_sum(tau_z);

      const chromosome &chrom = ligand_chromosomes[individual_index];
      const fp_type rad_beta  = deg_to_rad(chrom[4]);
      const fp_type rad_gamma = deg_to_rad(chrom[5]);

      const fp_type sin_beta  = std::sin(rad_beta);
      const fp_type cos_beta  = std::cos(rad_beta);
      const fp_type sin_gamma = std::sin(rad_gamma);
      const fp_type cos_gamma = std::cos(rad_gamma);

      const fp_type u_gamma[3] = {fp_type{0}, fp_type{0}, fp_type{1}};
      const fp_type u_beta[3]  = {-sin_gamma, cos_gamma, fp_type{0}};
      const fp_type u_alpha[3] = {cos_beta * cos_gamma, cos_beta * sin_gamma, -sin_beta};

      const fp_type g_rot_alpha = angle_scale * dot3(tau, u_alpha);
      const fp_type g_rot_beta  = angle_scale * dot3(tau, u_beta);
      const fp_type g_rot_gamma = angle_scale * dot3(tau, u_gamma);

      gradient &grad_l = gradients_l[individual_index];
      if (local_thread_id == 0) {
        grad_l[0] = g_trasl_x;
        grad_l[1] = g_trasl_y;
        grad_l[2] = g_trasl_z;
        grad_l[3] = g_rot_alpha;
        grad_l[4] = g_rot_beta;
        grad_l[5] = g_rot_gamma;
      }

      // ---------------------------------------------------------------------------------------------
      // 4. Torsion derivatives: every lane sums over its atoms, then warp all-reduce, lane 0 stores
      // ---------------------------------------------------------------------------------------------
      const fp_type gmean_x = g_trasl_x * inv_n;
      const fp_type gmean_y = g_trasl_y * inv_n;
      const fp_type gmean_z = g_trasl_z * inv_n;

      for (int t = 0; t < num_rotamers; ++t) {
        const int* frag_mask = fragments + t * num_atoms;

        const int a1 = frag_start_indices[t];
        const int a2 = frag_stop_indices[t];

        // axis is computed redundantly by every lane (same value, avoids a broadcast)
        fp_type axis[3];
        axis[0] = scratch_x_l[a2] - scratch_x_l[a1];
        axis[1] = scratch_y_l[a2] - scratch_y_l[a1];
        axis[2] = scratch_z_l[a2] - scratch_z_l[a1];

        const fp_type norm = std::sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
        const fp_type inv_norm = fp_type{1} / norm;
        axis[0] *= inv_norm;
        axis[1] *= inv_norm;
        axis[2] *= inv_norm;

        fp_type dE_dpsi{0};
        for (int i = local_thread_id; i < num_atoms; i += BLOCK_SIZE) {
          if (frag_mask[i] != 0) {
            const fp_type r[3] = {scratch_x_l[i] - scratch_x_l[a1],
                                  scratch_y_l[i] - scratch_y_l[a1],
                                  scratch_z_l[i] - scratch_z_l[a1]};
            fp_type temp[3];
            cross3(axis, r, temp);
            dE_dpsi += (dE_dx[i] - gmean_x) * temp[0] +
                       (dE_dy[i] - gmean_y) * temp[1] +
                       (dE_dz[i] - gmean_z) * temp[2];
          }
        }
        dE_dpsi = warp_allreduce_sum(dE_dpsi); // called by all lanes (no divergence before)
        if (local_thread_id == 0) grad_l[6 + t] = angle_scale * dE_dpsi;
      }

      // WAR hazard: the next individual overwrites dE_* in phase 1 while slower lanes could still be
      // reading them here
      __syncwarp();
    }
  }


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
                    (void*) &scores_b,
                    (void*) &converged_ligands_b};
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
    init_device(dev_id, minimum, maximum, center, map_index_xyz, grid_maps);
    
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
                    (void*) &ls_selected_b,
                    (void*) &converged_ligands_b};
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

  template<int MAX_ATOMS>
  batch_multiple get_evaluate_gradient_batch(const int device_id) {
    return get_kernel_batch_multiple_cuda<calc_gradient<MAX_ATOMS>>(device_id,
                                                                  BLOCK_SIZE,
                                                                  0,
                                                                  "adt_score::calc_gradient");
  }

  template<>
  batch_multiple get_adt_score_batch_multiple<queue_cuda>(const int atoms, std::shared_ptr<queue_cuda> q_b) {
    batch_multiple bucket_multiple{};
    const int device_id = q_b->get_id();
    constexpr_for<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>([&](const auto atom_index) {
      const auto n_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
      if (atoms == n_atoms)
        bucket_multiple = min_batch_multiple(get_evaluate_fitness_batch<n_atoms>(device_id), get_evaluate_gradient_batch<n_atoms>(device_id));
    });
    if (bucket_multiple.total_multiple() <= 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");
    return normalize_batch_multiple(bucket_multiple);
  }

} // namespace mudock
