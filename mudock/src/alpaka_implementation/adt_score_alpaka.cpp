#include <mudock/alpaka_implementation/adt_score_alpaka.hpp>

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

#ifndef MUDOCK_ALPAKA_BLOCK_SIZE
  #define MUDOCK_ALPAKA_BLOCK_SIZE 32
#endif

namespace mudock {
  namespace {
    ALPAKA_FN_ACC fp_type trilinear_interpolation_alpaka(const fp_type* map,
                                                         const fp_type* coeffs,
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
                                    const fp_type* scratch_x,
                                    const fp_type* scratch_y,
                                    const fp_type* scratch_z,
                                    const fp_type* vols_b,
                                    const fp_type* solpars_b,
                                    const fp_type* charges_b,
                                    const int* num_atoms_b,
                                    const int* num_rotamers_b,
                                    const int* num_nonbonds_b,
                                    const int* nonbond_a1_b,
                                    const int* nonbond_a2_b,
                                    const fp_type* nonbond_cA_b,
                                    const fp_type* nonbond_cB_b,
                                    const int* nonbond_xB_b,
                                    const fp_type* minimum,
                                    const fp_type* maximum,
                                    const fp_type* center,
                                    const int map_index_x,
                                    const int map_index_xy,
                                    const int map_index_xyz,
                                    const fp_type* grid_maps,
                                    const int* map_offsets_b,
                                    fp_type* scores_b) const {
        const int ligand_id = static_cast<int>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        const int thread_id = static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);

        // First implementation: keep the Alpaka backend correct and easy to test.
        // Parallelizing the atom/non-bond loops can follow once the backend is wired end-to-end.
        if (thread_id != 0) {
          return;
        }

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

        for (int scores_index = 0; scores_index < scores_per_ligand; ++scores_index) {
          const fp_type* ligand_x = l_scratch_x + scores_index * atom_stride;
          const fp_type* ligand_y = l_scratch_y + scores_index * atom_stride;
          const fp_type* ligand_z = l_scratch_z + scores_index * atom_stride;

          fp_type elect_total_trilinear = 0;
          fp_type emap_total_trilinear = 0;
          fp_type dmap_total_trilinear = 0;
          const fp_type* electro_map =
              grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
          const fp_type* desolv_map =
              grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

          for (int atom_index = 0; atom_index < MAX_ATOMS; ++atom_index) {
            if (atom_index < num_atoms) {
              fp_type coord[3]{ligand_x[atom_index], ligand_y[atom_index], ligand_z[atom_index]};

              if (coord[0] < minimum[0] || coord[0] > maximum[0] || coord[1] < minimum[1] ||
                  coord[1] > maximum[1] || coord[2] < minimum[2] || coord[2] > maximum[2]) {
                const auto diff_x = coord[0] - center[0];
                const auto diff_y = coord[1] - center[1];
                const auto diff_z = coord[2] - center[2];
                const fp_type distance_two = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
                const fp_type epenalty = distance_two * ENERGYPENALTY;
                elect_total_trilinear += epenalty;
                emap_total_trilinear += epenalty;
              } else {
                const auto atom_charge = l_charge[atom_index];
                const fp_type* atom_map = grid_maps + l_atom_map_offsets[atom_index];

                coord[0] = (coord[0] - minimum[0]) * inv_spacing;
                coord[1] = (coord[1] - minimum[1]) * inv_spacing;
                coord[2] = (coord[2] - minimum[2]) * inv_spacing;

                const int u0 = coord[0];
                const fp_type p0u = coord[0] - static_cast<fp_type>(u0);
                const fp_type p1u = fp_type{1} - p0u;

                const int v0 = coord[1];
                const fp_type p0v = coord[1] - static_cast<fp_type>(v0);
                const fp_type p1v = fp_type{1} - p0v;

                const int w0 = coord[2];
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
                    ::fabs(atom_charge);
              }
            }
          }

          fp_type elect_total_eintcal{0};
          fp_type emap_total_eintcal{0};
          fp_type dmap_total_eintcal{0};
          if (num_rotamers > 0) {
            for (int nonbond_index = 0; nonbond_index < num_nonbonds; ++nonbond_index) {
              const int a1 = l_nonbond_a1[nonbond_index];
              const int a2 = l_nonbond_a2[nonbond_index];

              const auto diff_x = ligand_x[a1] - ligand_x[a2];
              const auto diff_y = ligand_y[a1] - ligand_y[a2];
              const auto diff_z = ligand_z[a1] - ligand_z[a2];
              const fp_type distance_two = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
              const fp_type distance_two_clamp =
                  distance_two > RMIN_ELEC_SQUARE ? distance_two : RMIN_ELEC_SQUARE;
              const fp_type distance = ::sqrt(distance_two_clamp);

              const fp_type epsilon =
                  mehler_solmajer::A +
                  mehler_solmajer::B /
                      (fp_type{1} + mehler_solmajer::rk * ::exp(mehler_solmajer::lambda_B * distance));
              const fp_type r_dielectric = fp_type{1} / (distance * epsilon);
              const fp_type e_elec = l_charge[a1] * l_charge[a2] * ELECSCALE *
                                     autodock_parameters::coeff_estat * r_dielectric;
              elect_total_eintcal += e_elec;

              const fp_type nb_desolv =
                  (l_vol[a2] * (l_solpar[a1] + qsolpar * ::fabs(l_charge[a1])) +
                   l_vol[a1] * (l_solpar[a2] + qsolpar * ::fabs(l_charge[a2])));

              const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                       ::exp(fp_type{-0.5} / sigma_square * distance_two_clamp) *
                                       nb_desolv;
              dmap_total_eintcal += e_desolv;

              fp_type e_vdW_Hb{0};
              if (distance_two_clamp < nbc2) {
                const int xA = xA_default;
                const int xB = l_nonbond_xB[nonbond_index];

                if (xA != xB) {
                  const fp_type cA = l_nonbond_cA[nonbond_index];
                  const fp_type cB = l_nonbond_cB[nonbond_index];

                  const auto log_distance = ::log(distance);
                  const fp_type rA = ::exp(static_cast<fp_type>(xA) * log_distance);
                  const fp_type rB = ::exp(static_cast<fp_type>(xB) * log_distance);
                  const fp_type e = cA / rA - cB / rB;
                  e_vdW_Hb = EINTCLAMP < e ? EINTCLAMP : e;
                }
              }
              emap_total_eintcal += e_vdW_Hb;
            }
          }

          const fp_type tors_free_energy = num_rotamers * autodock_parameters::coeff_tors;
          scores_l[scores_index] = emap_total_eintcal + elect_total_eintcal + dmap_total_eintcal +
                                   emap_total_trilinear + elect_total_trilinear +
                                   dmap_total_trilinear + tors_free_energy;
        }
      }
    };
  } // namespace

  template<>
  void adt_score_kernel<queue_alpaka>::operator()() {
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
                                                   minimum,
                                                   maximum,
                                                   center,
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
  int get_adt_score_batch<queue_alpaka>(const int atoms,
                                        std::shared_ptr<queue_alpaka>,
                                        const size_t max_bucket_size) {
#ifdef MUDOCK_ADT_BUCKET_OVERRIDE
    const int capped = std::min<int>(MUDOCK_ADT_BUCKET_OVERRIDE, max_bucket_size);
    mudock::info("ALPAKA Bucket size for ",
                 atoms,
                 " atoms override -> ",
                 MUDOCK_ADT_BUCKET_OVERRIDE,
                 ", capped -> ",
                 capped);
    return capped;
#else
    mudock::info("ALPAKA Bucket size for ", atoms, " atoms -> ", max_bucket_size);
    return static_cast<int>(max_bucket_size);
#endif
  }
} // namespace mudock
