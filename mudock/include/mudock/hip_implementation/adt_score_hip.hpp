#pragma once

#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/adt_score_args.hpp>
#include <mudock/hip_implementation/queue_hip.hpp>
#include <mudock/utils.hpp>

#if defined(__HIPCC__) || defined(__CUDACC__)
  #include <mudock/hip_implementation/hip_utils.hpp>
#endif

namespace mudock {
#if defined(__HIPCC__) || defined(__CUDACC__)
  template<int MAX_ATOMS>
  __device__ __forceinline__ fp_type score_pose_hip(const adt_score_args& args) {
    const int lane          = threadIdx.x;
    const int ligand_id     = args.ligand_index;
    const int* __restrict__ num_atoms_b = args.num_atoms;
    const int* __restrict__ num_rotamers_b = args.num_rotamers;
    const int* __restrict__ num_nonbonds_b = args.num_nonbonds;
    const int num_atoms = num_atoms_b[ligand_id];
    const int num_rotamers = num_rotamers_b[ligand_id];
    const int nonbond_begin = num_nonbonds_b[ligand_id];
    const int num_nonbonds = num_nonbonds_b[ligand_id + 1] - nonbond_begin;
    const int atom_offset   = ligand_id * args.batch_atom_stride;
    const fp_type* __restrict__ x = args.x;
    const fp_type* __restrict__ y = args.y;
    const fp_type* __restrict__ z = args.z;
    const fp_type* __restrict__ vols = args.vols + atom_offset;
    const fp_type* __restrict__ solpars = args.solpars + atom_offset;
    const fp_type* __restrict__ charges = args.charges + atom_offset;
    const int* __restrict__ map_offsets = args.atom_map_offsets + atom_offset;
    const int* __restrict__ nonbond_a1 = args.nonbond_a1 + nonbond_begin;
    const int* __restrict__ nonbond_a2 = args.nonbond_a2 + nonbond_begin;
    const fp_type* __restrict__ nonbond_cA = args.nonbond_cA + nonbond_begin;
    const fp_type* __restrict__ nonbond_cB = args.nonbond_cB + nonbond_begin;
    const int* __restrict__ nonbond_xB = args.nonbond_xB + nonbond_begin;
    const fp_type* __restrict__ grid_maps = args.grid_maps;
    const fp_type* __restrict__ minimum = args.minimum;
    const fp_type* __restrict__ maximum = args.maximum;
    const fp_type* __restrict__ center = args.center;
    const fp_type* electro_map = grid_maps + args.map_index_xyz *
                                 static_cast<int>(autodock_grid_type::ELEC);
    const fp_type* desolv_map = grid_maps + args.map_index_xyz *
                                static_cast<int>(autodock_grid_type::DESOLV);
    const int map_size_y  = args.map_index_xy / args.map_index_x;
    const int map_size_z  = args.map_index_xyz / args.map_index_xy;
    const int last_cell_x = args.map_index_x - 2;
    const int last_cell_y = map_size_y - 2;
    const int last_cell_z = map_size_z - 2;

    fp_type electro_grid{0}, vdw_grid{0}, desolv_grid{0};
    MUDOCK_PRAGMA_TEMPLATE_UNROLL(MAX_ATOMS, BLOCK_SIZE)
    for (int atom = lane; atom < num_atoms; atom += BLOCK_SIZE) {
      const fp_type atom_x = x[atom], atom_y = y[atom], atom_z = z[atom];
      if (!isfinite(atom_x) || !isfinite(atom_y) || !isfinite(atom_z)) {
        electro_grid += EINTCLAMP;
        vdw_grid += EINTCLAMP;
      } else if (atom_x <= minimum[0] || atom_x >= maximum[0] || atom_y <= minimum[1] ||
                 atom_y >= maximum[1] || atom_z <= minimum[2] || atom_z >= maximum[2]) {
        const fp_type dx = atom_x - center[0], dy = atom_y - center[1], dz = atom_z - center[2];
        const fp_type penalty = (dx * dx + dy * dy + dz * dz) * ENERGYPENALTY;
        electro_grid += penalty;
        vdw_grid += penalty;
      } else {
        const fp_type gx = (atom_x - minimum[0]) * inv_spacing;
        const fp_type gy = (atom_y - minimum[1]) * inv_spacing;
        const fp_type gz = (atom_z - minimum[2]) * inv_spacing;
        const int ix = max(0, min(static_cast<int>(gx), last_cell_x));
        const int iy = max(0, min(static_cast<int>(gy), last_cell_y));
        const int iz = max(0, min(static_cast<int>(gz), last_cell_z));
        const fp_type fx = gx - static_cast<fp_type>(ix), fy = gy - static_cast<fp_type>(iy),
                      fz = gz - static_cast<fp_type>(iz);
        const fp_type px[2] = {fp_type{1} - fx, fx};
        const fp_type py[2] = {fp_type{1} - fy, fy};
        const fp_type pz[2] = {fp_type{1} - fz, fz};
        const fp_type coeffs[8] = {px[0] * py[0] * pz[0], px[0] * py[0] * pz[1],
                                   px[0] * py[1] * pz[0], px[0] * py[1] * pz[1],
                                   px[1] * py[0] * pz[0], px[1] * py[0] * pz[1],
                                   px[1] * py[1] * pz[0], px[1] * py[1] * pz[1]};
        const int base = args.map_index_xy * iz + args.map_index_x * iy + ix;
        const fp_type* atom_map = grid_maps + map_offsets[atom];
        fp_type elec{0}, desolv{0}, vdw{0};
        for (int corner = 0; corner < 8; ++corner) {
          const int offset = (corner & 4 ? 1 : 0) + (corner & 2 ? args.map_index_x : 0) +
                             (corner & 1 ? args.map_index_xy : 0);
          elec += coeffs[corner] * electro_map[base + offset];
          desolv += coeffs[corner] * desolv_map[base + offset];
          vdw += coeffs[corner] * atom_map[base + offset];
        }
        electro_grid += elec * charges[atom];
        desolv_grid += desolv * fabsf(charges[atom]);
        vdw_grid += vdw;
      }
    }

    fp_type electro_internal{0}, vdw_internal{0}, desolv_internal{0};
    if (num_rotamers > 0) {
      for (int pair = lane; pair < num_nonbonds; pair += BLOCK_SIZE) {
        const int a1 = nonbond_a1[pair], a2 = nonbond_a2[pair];
        const fp_type q1 = charges[a1], q2 = charges[a2];
        const fp_type dx = x[a1] - x[a2];
        const fp_type dy = y[a1] - y[a2];
        const fp_type dz = z[a1] - z[a2];
        const fp_type distance2 = dx * dx + dy * dy + dz * dz;
        const fp_type distance2_clamp = distance2 > RMIN_ELEC_SQUARE ? distance2 : RMIN_ELEC_SQUARE;
        const fp_type distance = sqrtf(distance2_clamp);
        const fp_type epsilon = mehler_solmajer::A + mehler_solmajer::B /
            (fp_type{1} + mehler_solmajer::rk * expf(mehler_solmajer::lambda_B * distance));
        electro_internal += q1 * q2 * ELECSCALE * autodock_parameters::coeff_estat /
                            (distance * epsilon);
        const fp_type nb_desolv = vols[a2] * (solpars[a1] + qsolpar * fabsf(q1)) +
                                  vols[a1] * (solpars[a2] + qsolpar * fabsf(q2));
        desolv_internal += autodock_parameters::coeff_desolv *
                           expf(fp_type{-0.5} / sigma_square * distance2_clamp) * nb_desolv;
        if (distance2_clamp < nbc2 && xA_default != nonbond_xB[pair]) {
          const fp_type log_distance = logf(distance);
          const fp_type rA = expf(static_cast<fp_type>(xA_default) * log_distance);
          const fp_type rB = expf(static_cast<fp_type>(nonbond_xB[pair]) * log_distance);
          fp_type energy = nonbond_cA[pair] / rA - nonbond_cB[pair] / rB;
          vdw_internal += energy < EINTCLAMP ? energy : EINTCLAMP;
        }
      }
    }

    fp_type energy = electro_grid + vdw_grid + desolv_grid + electro_internal + vdw_internal + desolv_internal;
    for (int offset = BLOCK_SIZE / 2; offset > 0; offset /= 2)
      energy += SHFL_DOWN(BITLANE_MASK, energy, offset, BLOCK_SIZE);
    if (lane == 0)
      energy += static_cast<fp_type>(num_rotamers) * autodock_parameters::coeff_tors;
    return SHFL(BITLANE_MASK, energy, 0, BLOCK_SIZE);
  }

  struct adt_score_device_kernel {
    using parameters = adt_score_args;

    template<int MAX_ATOMS>
    __device__ static fp_type score(const parameters& args) {
      return score_pose_hip<MAX_ATOMS>(args);
    }
  };
#endif

#if defined(__HIPCC__) || defined(__CUDACC__)
  template<int MAX_ATOMS>
  __global__ void calc_energy(const adt_score_args args);
#endif

  template<>
  void adt_score<queue_hip>::operator()();
} // namespace mudock
