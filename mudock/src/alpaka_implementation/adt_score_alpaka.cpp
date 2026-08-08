#include <mudock/alpaka_implementation/adt_score_alpaka.hpp>
#include <mudock/alpaka_implementation/invoke_kernel_alpaka.hpp>

#include <alpaka/alpaka.hpp>

#include <cmath>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>

#define FLATTENED_3D(x, y, z, index_x, index_xy) (index_xy * (z) + (y) * index_x + (x))



namespace mudock {
  
  ALPAKA_STATIC_ACC_MEM_CONSTANT alpaka::DevGlobal<TAcc, fp_type[3]> map_min_const;
  ALPAKA_STATIC_ACC_MEM_CONSTANT alpaka::DevGlobal<TAcc, fp_type[3]> map_max_const;
  ALPAKA_STATIC_ACC_MEM_CONSTANT alpaka::DevGlobal<TAcc, fp_type[3]> map_center_const;

  namespace {
    ALPAKA_FN_ACC ALPAKA_FN_INLINE fp_type trilinear_interpolation_alpaka(const fp_type* __restrict__ map,
                                                         const fp_type* __restrict__ coeffs,
                                                         const int map_index_x,
                                                         const int map_index_xy) {
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
    struct calc_energy {
      template<typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int atom_stride,
                                    const int scores_per_ligand,
                                    const fp_type* __restrict__ scratch_x,
                                    const fp_type* __restrict__ scratch_y,
                                    const fp_type* __restrict__ scratch_z,
                                    const fp_type* __restrict__ vols_b,
                                    const fp_type* __restrict__ solpars_b,
                                    const fp_type* __restrict__ charges_b,
                                    const int* num_atoms_b,
                                    const int* num_rotamers_b,
                                    const int* num_nonbonds_b,
                                    const int* __restrict__ nonbond_a1_b,
                                    const int* __restrict__ nonbond_a2_b,
                                    const fp_type* __restrict__ nonbond_cA_b,
                                    const fp_type* __restrict__ nonbond_cB_b,
                                    const int* __restrict__ nonbond_xB_b,
                                    const int map_index_x,
                                    const int map_index_xy,
                                    const int map_index_xyz,
                                    const fp_type* __restrict__ grid_maps,
                                    const int* __restrict__ map_offsets_b,
                                    fp_type* __restrict__ scores_b) const {
        const int ligand_id = static_cast<int>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        const int thread_id = static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);

        const int num_atoms = num_atoms_b[ligand_id];
        const int num_nonbonds = num_nonbonds_b[ligand_id + 1] - num_nonbonds_b[ligand_id];
        const int num_rotamers = num_rotamers_b[ligand_id];
        const int stride = ligand_id * atom_stride;

        const fp_type* l_scratch_x = scratch_x + stride * scores_per_ligand;
        const fp_type* l_scratch_y = scratch_y + stride * scores_per_ligand;
        const fp_type* l_scratch_z = scratch_z + stride * scores_per_ligand;
        const fp_type* l_vol = vols_b + stride;
        const fp_type* l_solpar = solpars_b + stride;
        const fp_type* l_charge = charges_b + stride;
        const int* l_atom_map_offsets = map_offsets_b + stride;
        const int* l_nonbond_a1 = nonbond_a1_b + num_nonbonds_b[ligand_id];
        const int* l_nonbond_a2 = nonbond_a2_b + num_nonbonds_b[ligand_id];
        const fp_type* l_nonbond_cA = nonbond_cA_b + num_nonbonds_b[ligand_id];
        const fp_type* l_nonbond_cB = nonbond_cB_b + num_nonbonds_b[ligand_id];
        const int* l_nonbond_xB = nonbond_xB_b + num_nonbonds_b[ligand_id];

        fp_type* scores_l = scores_b + ligand_id * scores_per_ligand;

        const fp_type l_minimum[3] = {map_min_const<TAcc>.get()[0], map_min_const<TAcc>.get()[1], map_min_const<TAcc>.get()[2]};
        const fp_type l_maximum[3] = {map_max_const<TAcc>.get()[0], map_max_const<TAcc>.get()[1], map_max_const<TAcc>.get()[2]};
        const fp_type l_center[3]  = {map_center_const<TAcc>.get()[0], map_center_const<TAcc>.get()[1], map_center_const<TAcc>.get()[2]};

        const fp_type* electro_map =
            grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
        const fp_type* desolv_map =
            grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

        for (int scores_index = 0; scores_index < scores_per_ligand; ++scores_index) {
          const fp_type* ligand_x = l_scratch_x + scores_index * atom_stride;
          const fp_type* ligand_y = l_scratch_y + scores_index * atom_stride;
          const fp_type* ligand_z = l_scratch_z + scores_index * atom_stride;

          fp_type elect_total_trilinear = 0;
          fp_type emap_total_trilinear = 0;
          fp_type dmap_total_trilinear = 0;

          ALPAKA_UNROLL(MUDOCK_ATOM_LOOP_UNROLL_FACTOR(MAX_ATOMS, MUDOCK_ALPAKA_BLOCK_SIZE))
          for (int atom_index = thread_id; atom_index < MAX_ATOMS; atom_index += MUDOCK_ALPAKA_BLOCK_SIZE) {
            if (atom_index < num_atoms) {
              fp_type coord[3]{ligand_x[atom_index], ligand_y[atom_index], ligand_z[atom_index]};

              if (coord[0] < l_minimum[0] || coord[0] > l_maximum[0] || coord[1] < l_minimum[1] ||
                  coord[1] > l_maximum[1] || coord[2] < l_minimum[2] || coord[2] > l_maximum[2]) {
                const auto diff_x = coord[0] - l_center[0];
                const auto diff_y = coord[1] - l_center[1];
                const auto diff_z = coord[2] - l_center[2];
                const fp_type distance_two = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
                const fp_type epenalty = distance_two * ENERGYPENALTY;
                elect_total_trilinear += epenalty;
                emap_total_trilinear += epenalty;
              } else {
                const auto atom_charge = l_charge[atom_index];
                const fp_type* atom_map = grid_maps + l_atom_map_offsets[atom_index];

                coord[0] = (coord[0] - l_minimum[0]) * inv_spacing;
                coord[1] = (coord[1] - l_minimum[1]) * inv_spacing;
                coord[2] = (coord[2] - l_minimum[2]) * inv_spacing;

                const int u0 = static_cast<int>(coord[0]);
                const fp_type p0u = coord[0] - static_cast<fp_type>(u0);
                const fp_type p1u = fp_type{1} - p0u;

                const int v0 = static_cast<int>(coord[1]);
                const fp_type p0v = coord[1] - static_cast<fp_type>(v0);
                const fp_type p1v = fp_type{1} - p0v;

                const int w0 = static_cast<int>(coord[2]);
                const fp_type p0w = coord[2] - static_cast<fp_type>(w0);
                const fp_type p1w = fp_type{1} - p0w;

                const fp_type pu[2] = {p1u, p0u};
                const fp_type pv[2] = {p1v, p0v};
                const fp_type pw[2] = {p1w, p0w};

                const fp_type coeffs[8] = {pu[0] * pv[0] * pw[0],
                                           pu[0] * pv[0] * pw[1],
                                           pu[0] * pv[1] * pw[0],
                                           pu[0] * pv[1] * pw[1],
                                           pu[1] * pv[0] * pw[0],
                                           pu[1] * pv[0] * pw[1],
                                           pu[1] * pv[1] * pw[0],
                                           pu[1] * pv[1] * pw[1]};

                const int base_index = FLATTENED_3D(u0, v0, w0, map_index_x, map_index_xy);
                elect_total_trilinear +=
                    trilinear_interpolation_alpaka(electro_map + base_index,
                                                   coeffs,
                                                   map_index_x,
                                                   map_index_xy) *
                    atom_charge;
                emap_total_trilinear +=
                    trilinear_interpolation_alpaka(atom_map + base_index,
                                                   coeffs,
                                                   map_index_x,
                                                   map_index_xy);
                dmap_total_trilinear +=
                    trilinear_interpolation_alpaka(desolv_map + base_index,
                                                   coeffs,
                                                   map_index_x,
                                                   map_index_xy) *
                    alpaka::math::abs(acc, atom_charge);
              }
            }
          }

          fp_type elect_total_eintcal{0};
          fp_type emap_total_eintcal{0};
          fp_type dmap_total_eintcal{0};
          if (num_rotamers > 0) {
            for (int nonbond_index = thread_id; nonbond_index < num_nonbonds; nonbond_index += MUDOCK_ALPAKA_BLOCK_SIZE) {
              const int a1 = l_nonbond_a1[nonbond_index];
              const int a2 = l_nonbond_a2[nonbond_index];

              const auto diff_x = ligand_x[a1] - ligand_x[a2];
              const auto diff_y = ligand_y[a1] - ligand_y[a2];
              const auto diff_z = ligand_z[a1] - ligand_z[a2];
              const fp_type distance_two = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
              const fp_type distance_two_clamp =
                  distance_two > RMIN_ELEC_SQUARE ? distance_two : RMIN_ELEC_SQUARE;
              const fp_type distance = alpaka::math::sqrt(acc, distance_two_clamp);

              const fp_type epsilon =
                  mehler_solmajer::A +
                  mehler_solmajer::B /
                      (fp_type{1} + mehler_solmajer::rk * alpaka::math::exp(acc, mehler_solmajer::lambda_B * distance));
              const fp_type r_dielectric = fp_type{1} / (distance * epsilon);
              const fp_type e_elec = l_charge[a1] * l_charge[a2] * ELECSCALE *
                                     autodock_parameters::coeff_estat * r_dielectric;
              elect_total_eintcal += e_elec;

              const fp_type nb_desolv =
                  (l_vol[a2] * (l_solpar[a1] + qsolpar * alpaka::math::abs(acc, l_charge[a1])) +
                   l_vol[a1] * (l_solpar[a2] + qsolpar * alpaka::math::abs(acc, l_charge[a2])));

              const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                       alpaka::math::exp(acc, fp_type{-0.5} / sigma_square * distance_two_clamp) *
                                       nb_desolv;
              dmap_total_eintcal += e_desolv;

              fp_type e_vdW_Hb{0};
              if (distance_two_clamp < nbc2) {
                const int xA = xA_default;
                const int xB = l_nonbond_xB[nonbond_index];

                if (xA != xB) {
                  const fp_type cA = l_nonbond_cA[nonbond_index];
                  const fp_type cB = l_nonbond_cB[nonbond_index];

                  const auto log_distance = alpaka::math::log(acc, distance);
                  const fp_type rA = alpaka::math::exp(acc, static_cast<fp_type>(xA) * log_distance);
                  const fp_type rB = alpaka::math::exp(acc, static_cast<fp_type>(xB) * log_distance);
                  const fp_type e = cA / rA - cB / rB;
                  e_vdW_Hb = EINTCLAMP < e ? EINTCLAMP : e;
                }
              }
              emap_total_eintcal += e_vdW_Hb;
            }
          }

          fp_type total_energy = emap_total_eintcal + elect_total_eintcal + dmap_total_eintcal +
                                 emap_total_trilinear + elect_total_trilinear +
                                 dmap_total_trilinear;

          ALPAKA_UNROLL(MUDOCK_UNROLL_FACTOR)
          for (int offset = MUDOCK_ALPAKA_BLOCK_SIZE / 2; offset > 0; offset /= 2) {
            total_energy += alpaka::warp::shfl_down(acc, total_energy, offset, MUDOCK_ALPAKA_BLOCK_SIZE);
          }

          if (thread_id == 0) {
            const fp_type tors_free_energy =
                static_cast<fp_type>(num_rotamers) * autodock_parameters::coeff_tors;
            scores_l[scores_index] = total_energy + tors_free_energy;
          }
        }
      }
    };
  } // namespace

  template<>
  void adt_score_kernel<queue_alpaka>::operator()() {
    const auto extent = alpaka::Vec<alpaka_backend::dim, alpaka_backend::idx>{3u};
    auto view_min = alpaka::createView(q->native_device(), minimum, extent);
    auto view_max = alpaka::createView(q->native_device(), maximum, extent);
    auto view_center = alpaka::createView(q->native_device(), center, extent);

    alpaka::memcpy(q->native_queue(), map_min_const<alpaka_backend::acc>, view_min, extent);
    alpaka::memcpy(q->native_queue(), map_max_const<alpaka_backend::acc>, view_max, extent);
    alpaka::memcpy(q->native_queue(), map_center_const<alpaka_backend::acc>, view_center, extent);

    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->invoke_kernel<calc_energy<max_atoms>>(batch_ligands,
                                                   MUDOCK_ALPAKA_BLOCK_SIZE,
                                                   batch_atoms,
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
                                                   map_index_x,
                                                   map_index_xy,
                                                   map_index_xyz,
                                                   grid_maps,
                                                   map_offsets_b,
                                                   scores_b);
        },
        batch_atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

  template<>
  batch_multiple get_adt_score_batch_multiple<queue_alpaka>(const int atoms,
                                                            std::shared_ptr<queue_alpaka> q_b) {
    const auto& dev = q_b->native_device();
    const int num_sms = static_cast<int>(alpaka::getAccDevProps<alpaka_backend::acc>(dev).m_multiProcessorCount);

#if defined(MUDOCK_ALPAKA_BACKEND_SERIAL) || defined(MUDOCK_ALPAKA_BACKEND_TBB) || defined(MUDOCK_ALPAKA_BACKEND_OMP2)
    const int blocks_per_sm = 1;
#else
    const int blocks_per_sm = 16;
#endif

    mudock::info("ALPAKA ADT batch multiple for ", atoms, " atoms -> ", blocks_per_sm * num_sms);
    return {blocks_per_sm, num_sms};
  }
} // namespace mudock
