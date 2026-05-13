#include "mudock/compute/queue.hpp"

#include <memory>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/compute/adt_score_kernel.hpp>
#include <mudock/compute/devices_memory.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/hip_implementation/adt_score_hip.hpp>
#include <mudock/hip_implementation/hip_utils.hpp>
#include <mudock/log.hpp>
#include <mudock/hip_implementation/hip_texture.hpp>
#include <mudock/hip_implementation/hip_utils.hpp>
#include <mudock/hip_implementation/queue_hip.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>

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

  device_memory_array<k_max_devices, hip_texture_devices> hip_texture_memory;
  device_memory_array<k_max_devices, fp_type> hip_constant_memory;

  void init_device(const int dev,
                   const fp_type* map_min,
                   const fp_type* map_max,
                   const fp_type* map_center,
                   const int map_index_x,
                   const int map_index_xy,
                   const int map_index_xyz,
                   const fp_type* map_grids) {
    // Thread-safe, exactly-once init per device:
    hip_texture_memory.init(dev, map_index_xyz, num_autodock_grids(), map_grids);
    hip_constant_memory.init(
        dev,
        std::function<std::unique_ptr<fp_type>()>([&]() {
          MUDOCK_CHECK(hipSetDevice(dev));
          MUDOCK_CHECK(
              hipMemcpyToSymbol(map_min_const, map_min, 3 * sizeof(fp_type), 0, hipMemcpyDeviceToDevice));
          MUDOCK_CHECK(
              hipMemcpyToSymbol(map_max_const, map_max, 3 * sizeof(fp_type), 0, hipMemcpyDeviceToDevice));
          MUDOCK_CHECK(hipMemcpyToSymbol(map_center_const,
                                         map_center,
                                         3 * sizeof(fp_type),
                                         0,
                                         hipMemcpyDeviceToDevice));
          MUDOCK_CHECK(hipDeviceSynchronize());
          // TODO implement RAI to release constant memory
          return std::make_unique<fp_type>(0);
        }));
  }

  __device__ __forceinline__ fp_type trilinear_interpolation_hip(const fp_type* __restrict__ map,
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
      for (int i = 0; i < MAX_ATOMS; i += BLOCK_SIZE) {
        const int atom_index = i + threadIdx.x;
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
            coord_tex[0]            = (coord_tex[0] - map_min_const[0]) * inv_spacing,
            coord_tex[1]            = (coord_tex[1] - map_min_const[1]) * inv_spacing;
            coord_tex[2]            = (coord_tex[2] - map_min_const[2]) * inv_spacing;
            const auto& charge      = l_charge[atom_index];
            const fp_type* atom_map = grid_maps + l_atom_tex_indexes[atom_index];

            //  TODO check approximations with in hardware interpolation
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
            elect_total_trilinear +=
                trilinear_interpolation_hip(electro_map + base_index, coeffs, map_index_x, map_index_xy) *
                charge;
            dmap_total_trilinear +=
                trilinear_interpolation_hip(desolv_map + base_index, coeffs, map_index_x, map_index_xy) *
                fabsf(charge);
            emap_total_trilinear +=
                trilinear_interpolation_hip(atom_map + base_index, coeffs, map_index_x, map_index_xy);
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
        total_energy += SHFL_DOWN(BITLANE_MASK, total_energy, offset, BLOCK_SIZE);
      }

      if (local_thread_id == 0) {
        const fp_type tors_free_energy = num_rotamers * autodock_parameters::coeff_tors;
        scores_l[scores_index]         = total_energy + tors_free_energy;
      }
    }
  }

  template<>
  void adt_score_kernel<queue_hip>::operator()() {
    const int dev_id = q->get_id();
    init_device(dev_id, minimum, maximum, center, map_index_x, map_index_xy, map_index_xyz, grid_maps);

    void* args[] = {(void*) &batch_atoms,    (void*) &scores_per_ligand,
                    (void*) &x_scratch_b,    (void*) &y_scratch_b,
                    (void*) &z_scratch_b,    (void*) &vols_b,
                    (void*) &solpars_b,      (void*) &charges_b,
                    (void*) &num_atoms_b,    (void*) &num_rotamers_b,
                    (void*) &num_nonbonds_b, (void*) &nonbond_a1_b,
                    (void*) &nonbond_a2_b,   (void*) &nonbond_cA_b,
                    (void*) &nonbond_cB_b,   (void*) &nonbond_xB_b,
                    (void*) &map_index_x,    (void*) &map_index_xy,
                    (void*) &map_index_xyz,  (void*) &hip_texture_memory.v[dev_id].data->tex_dev,
                    (void*) &map_offsets_b,  (void*) &scores_b};
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->launch_kernel((void*) calc_energy<max_atoms>, args, batch_ligands, BLOCK_SIZE);
        },
        batch_atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());

  }; // namespace mudock

  template<int MAX_ATOMS>
  batch_multiple get_evaluate_fitness_batch(const int device_id) {
    return get_kernel_batch_multiple_hip<calc_energy<MAX_ATOMS>>(device_id,
                                                                 BLOCK_SIZE,
                                                                 0,
                                                                 "adt_score::calc_energy");
  }

  template<>
  batch_multiple get_adt_score_batch_multiple<queue_hip>(const int atoms, std::shared_ptr<queue_hip> q_b) {
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
