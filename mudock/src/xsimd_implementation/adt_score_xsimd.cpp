#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/xsimd_implementation/adt_score_xsimd.hpp>
#include <xsimd/xsimd.hpp>

namespace mudock {
  namespace {
    template<typename T, typename V, typename VI>
    V trilinear_interpolation_vectorized_xsimd(const T* __restrict__ map,
                                               const VI& base_index,
                                               const V& p1u,
                                               const V& p1v,
                                               const V& p1w,
                                               const V& p0u,
                                               const V& p0v,
                                               const V& p0w,
                                               const VI& map_index_plus_one_vec,
                                               const VI& map_index_x_vec,
                                               const VI& map_index_x_plus_one_vec,
                                               const VI& map_index_xy_vec,
                                               const VI& map_index_xy_plus_one_vec,
                                               const VI& map_index_x_xy_vec,
                                               const VI& map_index_x_xy_plus_one_vec) {
    using batch_type = V;

    // Precompute products
    const auto p1v_p1w = p1v * p1w;
    const auto p1v_p0w = p1v * p0w;
    const auto p0v_p1w = p0v * p1w;
    const auto p0v_p0w = p0v * p0w;

    // Initialize value
    const auto zeros = batch_type(0);
    auto value       = zeros;

    // Gather and accumulate values
    value = xsimd::fma(p1u * p1v_p1w, batch_type::gather(map, base_index), value);
    value = xsimd::fma(p0u * p1v_p1w, batch_type::gather(map, base_index + map_index_plus_one_vec), value);
    value = xsimd::fma(p1u * p1v_p0w, batch_type::gather(map, base_index + map_index_xy_vec), value);
    value = xsimd::fma(p0u * p1v_p0w, batch_type::gather(map, base_index + map_index_xy_plus_one_vec), value);
    value = xsimd::fma(p1u * p0v_p1w, batch_type::gather(map, base_index + map_index_x_vec), value);
    value = xsimd::fma(p0u * p0v_p1w, batch_type::gather(map, base_index + map_index_x_plus_one_vec), value);
    value = xsimd::fma(p1u * p0v_p0w, batch_type::gather(map, base_index + map_index_x_xy_vec), value);
    value =
        xsimd::fma(p0u * p0v_p0w, batch_type::gather(map, base_index + map_index_x_xy_plus_one_vec), value);

    return value;
    }
    void calc_energy(const int batch_atoms,
                     const int batch_ligands,
                     const int scores_per_ligand,
                     const fp_type* __restrict__ x_scratch_b,
                     const fp_type* __restrict__ y_scratch_b,
                     const fp_type* __restrict__ z_scratch_b,
                     const fp_type* __restrict__ vols_b,
                     const fp_type* __restrict__ solpars_b,
                     const fp_type* __restrict__ charges_b,
                     const int* __restrict__ num_atoms_b,
                     const int* __restrict__ num_rotamers_b,
                     const int* __restrict__ num_nonbonds_b,
                     const int* __restrict__ nonbond_a1_b,
                     const int* __restrict__ nonbond_a2_b,
                     const fp_type* __restrict__ nonbond_cA_b,
                     const fp_type* __restrict__ nonbond_cB_b,
                     const int* __restrict__ nonbond_xB_b,
                     const fp_type* __restrict__ grid_maps,
                     const fp_type* __restrict__ minimum,
                     const fp_type* __restrict__ maximum,
                     const fp_type* __restrict__ center,
                     const int* __restrict__ map_offsets_b,
                     const int map_index_x,
                     const int map_index_xy,
                     const int map_index_xyz,
                     fp_type* __restrict__ scores_b) {
    using batch_type = xsimd::batch<fp_type>;
    using batch_int  = xsimd::batch<int>;
    using mask_type  = typename batch_type::batch_bool_type;

    constexpr std::size_t simd_size     = batch_type::size;
    constexpr std::size_t simd_size_int = batch_int::size;
    constexpr auto full_bitmask         = std::numeric_limits<uint>::max() >> (sizeof(uint) * 8 - simd_size);
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      const int atom_stride  = ligand_index * batch_atoms;
      const int num_atoms    = num_atoms_b[ligand_index];
      const int num_nonbonds = num_nonbonds_b[ligand_index + 1] - num_nonbonds_b[ligand_index];
      const int num_rotamers = num_rotamers_b[ligand_index];

      const fp_type* __restrict__ scratch_x = x_scratch_b + ligand_index * atom_stride * scores_per_ligand;
      const fp_type* __restrict__ scratch_y = y_scratch_b + ligand_index * atom_stride * scores_per_ligand;
      const fp_type* __restrict__ scratch_z = z_scratch_b + ligand_index * atom_stride * scores_per_ligand;
      const fp_type* __restrict__ vol_l     = vols_b + ligand_index * atom_stride;
      const fp_type* __restrict__ solpar_l  = solpars_b + ligand_index * atom_stride;
      const fp_type* __restrict__ charge_l  = charges_b + ligand_index * atom_stride;
      const int* __restrict__ map_offsets_l = map_offsets_b + ligand_index * atom_stride;
      const int* __restrict__ nonbond_a1_l  = nonbond_a1_b + num_nonbonds_b[ligand_index];
      const int* __restrict__ nonbond_a2_l  = nonbond_a2_b + num_nonbonds_b[ligand_index];
      const fp_type* nonbond_cA_l           = nonbond_cA_b + num_nonbonds_b[ligand_index];
      const fp_type* nonbond_cB_l           = nonbond_cB_b + num_nonbonds_b[ligand_index];
      const int* nonbond_xB_l               = nonbond_xB_b + num_nonbonds_b[ligand_index];

      fp_type* __restrict__ scores_l = scores_b + ligand_index * scores_per_ligand;
      for (int scores_index = 0; scores_index < scores_per_ligand; ++scores_index) {
        const fp_type* __restrict__ scratch_x_l = scratch_x + scores_index * batch_atoms;
        const fp_type* __restrict__ scratch_y_l = scratch_y + scores_index * batch_atoms;
        const fp_type* __restrict__ scratch_z_l = scratch_z + scores_index * batch_atoms;

        const auto num_atom_loops = static_cast<size_t>((num_atoms + simd_size - 1) / simd_size);
        const auto num_non_bond_loops =
            static_cast<size_t>((num_nonbonds + simd_size_int - 1) / simd_size_int);
        const fp_type* electro_map = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
        const fp_type* desolv_map  = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

        fp_type elect_total_trilinear = 0;
        fp_type emap_total_trilinear  = 0;
        fp_type dmap_total_trilinear  = 0;

        // Set up constants for SIMD operations
        const auto min_x = batch_type(minimum[0]);
        const auto max_x = batch_type(maximum[0]);
        const auto min_y = batch_type(minimum[1]);
        const auto max_y = batch_type(maximum[1]);
        const auto min_z = batch_type(minimum[2]);
        const auto max_z = batch_type(maximum[2]);

        const auto center_x = batch_type(center[0]);
        const auto center_y = batch_type(center[1]);
        const auto center_z = batch_type(center[2]);

        const auto inv_spacing_vec             = batch_type(inv_spacing);
        const auto one                         = batch_type(fp_type(1.0));
        const auto penalty_vec                 = batch_type(ENERGYPENALTY);
        const auto default_offset              = batch_int(0);
        const auto map_index_plus_one_vec      = batch_int(1);
        const auto map_index_x_vec             = batch_int(map_index_x);
        const auto map_index_x_plus_one_vec    = batch_int(map_index_x + 1);
        const auto map_index_xy_vec            = batch_int(map_index_xy);
        const auto map_index_xy_plus_one_vec   = batch_int(map_index_xy + 1);
        const auto map_index_x_xy_vec          = batch_int(map_index_x + map_index_xy);
        const auto map_index_x_xy_plus_one_vec = batch_int(map_index_x + map_index_xy + 1);

        for (std::size_t index = 0; index < num_atom_loops * simd_size; index += simd_size) {
          const auto remaining =
              static_cast<int>(index + simd_size) < num_atoms ? simd_size : (num_atoms - index) % simd_size;
          const auto remaining_mask =
              xsimd::batch_bool<fp_type>::from_mask(full_bitmask >> (simd_size - remaining));

          // Load coordinates
          auto x = xsimd::load_unaligned(scratch_x_l + index);
          auto y = xsimd::select(remaining_mask, xsimd::load_unaligned(scratch_y_l + index), batch_type{0});
          auto z = xsimd::load_unaligned(scratch_z_l + index);
          const auto atom_charge = xsimd::load_unaligned(charge_l + index);

          // Bounds check
          const mask_type outside_x = (x < min_x) | (x > max_x);
          const mask_type outside_y = (y < min_y) | (x > max_y);
          const mask_type outside_z = (z < min_z) | (z > max_z);
          const mask_type outside   = (outside_x | outside_y | outside_z) & remaining_mask;

          // Handle atoms outside boundaries
          if (xsimd::any(outside)) {
            const auto dx      = x - center_x;
            const auto dy      = y - center_y;
            const auto dz      = z - center_z;
            const auto dist    = dx * dx + dy * dy + dz * dz;
            const auto penalty = dist * penalty_vec;

            elect_total_trilinear += xsimd::reduce_add(xsimd::select(outside, penalty, batch_type(0)));
            emap_total_trilinear += xsimd::reduce_add(xsimd::select(outside, penalty, batch_type(0)));
          }

          // Handle atoms inside boundaries
          const mask_type inside = ~outside & remaining_mask;

          if (xsimd::any(inside)) {
            // Calculate interpolation coordinates
            x = (x - min_x) * inv_spacing_vec;
            y = (y - min_y) * inv_spacing_vec;
            z = (z - min_z) * inv_spacing_vec;

            // Load map offsets
            auto atom_map_offsets = xsimd::load_unaligned(map_offsets_l + index);

            // Decompose coordinates
            auto u0 = xsimd::to_int(x);
            auto v0 = xsimd::to_int(y);
            auto w0 = xsimd::to_int(z);

            auto p0u = x - xsimd::to_float(u0);
            auto p0v = y - xsimd::to_float(v0);
            auto p0w = z - xsimd::to_float(w0);

            auto p1u = one - p0u;
            auto p1v = one - p0v;
            auto p1w = one - p0w;

            // Compute base indices
            auto base_default_index = v0 * map_index_x_vec + w0 * map_index_xy_vec + u0 + default_offset;

            // Compute interpolated values
            elect_total_trilinear += xsimd::reduce_add(
                xsimd::select(inside,
                              trilinear_interpolation_vectorized_xsimd(electro_map,
                                                                       base_default_index,
                                                                       p1u,
                                                                       p1v,
                                                                       p1w,
                                                                       p0u,
                                                                       p0v,
                                                                       p0w,
                                                                       map_index_plus_one_vec,
                                                                       map_index_x_vec,
                                                                       map_index_x_plus_one_vec,
                                                                       map_index_xy_vec,
                                                                       map_index_xy_plus_one_vec,
                                                                       map_index_x_xy_vec,
                                                                       map_index_x_xy_plus_one_vec) *
                                  atom_charge,
                              batch_type(0)));

            dmap_total_trilinear += xsimd::reduce_add(
                xsimd::select(inside,
                              trilinear_interpolation_vectorized_xsimd(desolv_map,
                                                                       base_default_index,
                                                                       p1u,
                                                                       p1v,
                                                                       p1w,
                                                                       p0u,
                                                                       p0v,
                                                                       p0w,
                                                                       map_index_plus_one_vec,
                                                                       map_index_x_vec,
                                                                       map_index_x_plus_one_vec,
                                                                       map_index_xy_vec,
                                                                       map_index_xy_plus_one_vec,
                                                                       map_index_x_xy_vec,
                                                                       map_index_x_xy_plus_one_vec) *
                                  xsimd::abs(atom_charge),
                              batch_type(0)));

            auto base_atom_index = v0 * map_index_x_vec + w0 * map_index_xy_vec + u0 + atom_map_offsets;

            emap_total_trilinear += xsimd::reduce_add(
                xsimd::select(inside,
                              trilinear_interpolation_vectorized_xsimd(grid_maps,
                                                                       base_atom_index,
                                                                       p1u,
                                                                       p1v,
                                                                       p1w,
                                                                       p0u,
                                                                       p0v,
                                                                       p0w,
                                                                       map_index_plus_one_vec,
                                                                       map_index_x_vec,
                                                                       map_index_x_plus_one_vec,
                                                                       map_index_xy_vec,
                                                                       map_index_xy_plus_one_vec,
                                                                       map_index_x_xy_vec,
                                                                       map_index_x_xy_plus_one_vec),
                              batch_type(0)));
          }
        }

        fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
        if (num_rotamers > 0) {
          // Constants setup
          const auto ms_A_vec                    = batch_type(mehler_solmajer::A);
          const auto ms_B_vec                    = batch_type(mehler_solmajer::B);
          const auto ms_rk_vec                   = batch_type(mehler_solmajer::rk);
          const auto ms_lambda_B_vec             = batch_type(mehler_solmajer::lambda_B);
          const auto elec_scale                  = batch_type(ELECSCALE);
          const auto coeff_estat_vec             = batch_type(autodock_parameters::coeff_estat);
          const auto nbc2_vec                    = batch_type(nbc2);
          const auto half_fp_vec                 = batch_type(-0.5);
          const auto xA_vec                      = batch_type(12);
          const auto eintclamp_vec               = batch_type(EINTCLAMP);
          const auto rmin_elec_square_vec        = batch_type(RMIN_ELEC_SQUARE);
          const auto qsolpar_vec                 = batch_type(0.01097);
          const auto coeff_desolv_vec            = batch_type(autodock_parameters::coeff_desolv);
          const auto reciprocal_sigma_square_vec = batch_type(1.0) / batch_type(sigma_square);

          for (size_t i = 0; i < num_non_bond_loops * simd_size; i += simd_size) {
            const auto remaining =
                static_cast<int>(i + simd_size) < num_nonbonds ? simd_size : (num_nonbonds - i) % simd_size;
            const auto remaining_mask =
                xsimd::batch_bool<batch_type::value_type>::from_mask(full_bitmask >> (simd_size - remaining));

            const auto a1 = xsimd::load_unaligned(nonbond_a1_l + i);
            const auto a2 = xsimd::load_unaligned(nonbond_a2_l + i);

            const auto diff_x = batch_type::gather(scratch_x_l, a1) - batch_type::gather(scratch_x_l, a2);
            const auto diff_y = batch_type::gather(scratch_y_l, a1) - batch_type::gather(scratch_y_l, a2);
            const auto diff_z = batch_type::gather(scratch_z_l, a1) - batch_type::gather(scratch_z_l, a2);

            // Calculate squared distances
            const auto distance_two         = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
            const auto clamped_distance_two = xsimd::max(rmin_elec_square_vec, distance_two);
            const auto distance             = clamped_distance_two * xsimd::rsqrt(clamped_distance_two);

            // Calculate Electrostatic Energy
            const auto epsilon =
                ms_B_vec * (batch_type(1.0) /
                            (batch_type(1.0) + ms_rk_vec * xsimd::exp(ms_lambda_B_vec * distance))) +
                ms_A_vec;

            const auto r_dielectric = batch_type(1.0) / (distance * epsilon);

            const auto charge_a1 = batch_type::gather(charge_l, a1);
            const auto charge_a2 = batch_type::gather(charge_l, a2);

            const auto e_elec = charge_a1 * charge_a2 * elec_scale * coeff_estat_vec * r_dielectric;

            elect_total_eintcal += xsimd::reduce_add(xsimd::select(remaining_mask, e_elec, batch_type{0}));

            // Calculate desolv
            const auto ligand_vol_a1_v    = batch_type::gather(vol_l, a1);
            const auto ligand_vol_a2_v    = batch_type::gather(vol_l, a2);
            const auto ligand_solpar_a1_v = batch_type::gather(solpar_l, a1);
            const auto ligand_solpar_a2_v = batch_type::gather(solpar_l, a2);

            const auto nb_desolv =
                ligand_vol_a2_v * (ligand_solpar_a1_v + qsolpar_vec * xsimd::abs(charge_a1)) +
                ligand_vol_a1_v * (ligand_solpar_a2_v + qsolpar_vec * xsimd::abs(charge_a2));

            const auto e_desolv =
                coeff_desolv_vec *
                xsimd::exp(half_fp_vec * reciprocal_sigma_square_vec * clamped_distance_two) * nb_desolv;

            dmap_total_eintcal += xsimd::reduce_add(xsimd::select(remaining_mask, e_desolv, batch_type{0}));

            // Calculate vdW/Hb
            auto low_distance_mask = (clamped_distance_two < nbc2_vec);
            batch_type e_vdW_Hb(0.0);
            if (xsimd::any(low_distance_mask)) {
              const auto xB_vec   = batch_type::load_unaligned(nonbond_xB_l + i);
              const auto xab_cond = (xA_vec != xB_vec);

              if (xsimd::any(xab_cond)) {
                const auto cA = batch_type::load_unaligned(nonbond_cA_l + i);
                const auto cB = batch_type::load_unaligned(nonbond_cB_l + i);

                const auto log_distance = xsimd::log(distance);
                const auto rA           = xsimd::exp(batch_type(xA_vec) * log_distance);
                const auto rB           = xsimd::exp(batch_type(xB_vec) * log_distance);

                e_vdW_Hb = xsimd::min(eintclamp_vec, cA / rA - cB / rB);

                // Apply xab_cond mask
                e_vdW_Hb = xsimd::select(xab_cond, e_vdW_Hb, batch_type(0.0));
              }
            }

            // Accumulate all lanes
            emap_total_eintcal += xsimd::reduce_add(xsimd::select(remaining_mask, e_vdW_Hb, batch_type{0}));
          }
        }

        const fp_type tors_free_energy = num_rotamers * autodock_parameters::coeff_tors;
        const fp_type total_trilinear  = emap_total_trilinear + elect_total_trilinear + dmap_total_trilinear;
        const fp_type total_eintcal    = emap_total_eintcal + elect_total_eintcal + dmap_total_eintcal;

        scores_l[scores_index] = total_trilinear + total_eintcal + tors_free_energy;
      }
    }
    }
  } // namespace
  template<>
  void adt_score_kernel<queue_xsimd>::operator()() {
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
