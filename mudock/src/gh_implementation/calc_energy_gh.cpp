#include <hwy/contrib/math/math-inl.h>
#include <hwy/highway.h>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/cpp_implementation/calc_energy_gh.hpp>
#include <mudock/type_alias.hpp>
#include <stdint.h>

namespace mudock {
  template<typename VM, typename V, typename VI, typename T>
  inline V trilinear_interpolation_vectorized(const T* __restrict__ map,
                                              const VI& base_index,
                                              const VM inside,
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
    const HWY_FULL(T) d;

    // Precompute flattened indices
    auto value = Zero(d);

    const auto p1v_p1w = Mul(p1v, p1w);
    const auto p1v_p0w = Mul(p1v, p0w);
    const auto p0v_p1w = Mul(p0v, p1w);
    const auto p0v_p0w = Mul(p0v, p0w);

    // Gather and accumulate the values based on offsets
    value = MulAdd(Mul(p1u, p1v_p1w), MaskedGatherIndex(inside, d, map, base_index), value);
    value = MulAdd(Mul(p0u, p1v_p1w),
                   MaskedGatherIndex(inside, d, map, Add(base_index, map_index_plus_one_vec)),
                   value);
    value = MulAdd(Mul(p1u, p1v_p0w),
                   MaskedGatherIndex(inside, d, map, Add(base_index, map_index_xy_vec)),
                   value);
    value = MulAdd(Mul(p0u, p1v_p0w),
                   MaskedGatherIndex(inside, d, map, Add(base_index, map_index_xy_plus_one_vec)),
                   value);
    value =
        MulAdd(Mul(p1u, p0v_p1w), MaskedGatherIndex(inside, d, map, Add(base_index, map_index_x_vec)), value);
    value = MulAdd(Mul(p0u, p0v_p1w),
                   MaskedGatherIndex(inside, d, map, Add(base_index, map_index_x_plus_one_vec)),
                   value);
    value = MulAdd(Mul(p1u, p0v_p0w),
                   MaskedGatherIndex(inside, d, map, Add(base_index, map_index_x_xy_vec)),
                   value);
    value = MulAdd(Mul(p0u, p0v_p0w),
                   MaskedGatherIndex(inside, d, map, Add(base_index, map_index_x_xy_plus_one_vec)),
                   value);

    return value;
  }

  template<>
  fp_type calc_energy<cpu_vectorization::GH>(const fp_type* __restrict__ ligand_x,
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
                                             const int map_index_xyz,
                                             const fp_type* __restrict__ grid_maps) {
    const HWY_FULL(fp_type) d;
    const HWY_FULL(int) di;
    const auto num_atom_loops     = static_cast<size_t>((num_atoms + Lanes(d) - 1) / Lanes(d));
    const auto num_non_bond_loops = static_cast<size_t>((num_nonbond + Lanes(di) - 1) / Lanes(di));
    const fp_type* electro_map    = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
    const fp_type* desolv_map     = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

    static_assert(Lanes(d) == Lanes(di));

    fp_type elect_total_trilinear = 0;
    fp_type emap_total_trilinear  = 0;
    fp_type dmap_total_trilinear  = 0;

    // Set up constants for SIMD operations
    const auto min_x = Set(d, minimum[0]);
    const auto max_x = Set(d, maximum[0]);
    const auto min_y = Set(d, minimum[1]);
    const auto max_y = Set(d, maximum[1]);
    const auto min_z = Set(d, minimum[2]);
    const auto max_z = Set(d, maximum[2]);

    const auto center_x = Set(d, center[0]);
    const auto center_y = Set(d, center[1]);
    const auto center_z = Set(d, center[2]);

    const auto inv_spacing_vec             = Set(d, inv_spacing);
    const auto one                         = Set(d, fp_type{1.0});
    const auto penalty_vec                 = Set(d, ENERGYPENALTY);
    const auto default_offset              = Set(di, 0);
    const auto map_index_plus_one_vec      = Set(di, 1);
    const auto map_index_x_vec             = Set(di, map_index_x);
    const auto map_index_x_plus_one_vec    = Set(di, map_index_x + 1);
    const auto map_index_xy_vec            = Set(di, map_index_xy);
    const auto map_index_xy_plus_one_vec   = Set(di, map_index_xy + 1);
    const auto map_index_x_xy_vec          = Set(di, map_index_x + map_index_xy);
    const auto map_index_x_xy_plus_one_vec = Set(di, map_index_x + map_index_xy + 1);

    for (size_t index = 0; index < num_atom_loops * Lanes(di); index += Lanes(di)) {
      const auto remaining = num_atoms - index;
      // Load the x, y, z coordinates in a SIMD fashion
      auto x                 = LoadN(d, ligand_x + index, remaining);
      auto y                 = LoadN(d, ligand_y + index, remaining);
      auto z                 = LoadN(d, ligand_z + index, remaining);
      const auto valid_coord = FirstN(d, remaining);
      const auto atom_charge = LoadN(d, ligand_charge + index, remaining);

      // Bounds check for each atom in vectorized form
      const auto outside_x = Or(Lt(x, min_x), Gt(x, max_x));
      const auto outside_y = Or(Lt(y, min_y), Gt(y, max_y));
      const auto outside_z = Or(Lt(z, min_z), Gt(z, max_z));
      const auto outside   = And(Or(outside_x, Or(outside_y, outside_z)), valid_coord);

      // For atoms outside the boundaries
      if (FindFirstTrue(d, outside) != -1) {
        const auto dx      = Sub(x, center_x);
        const auto dy      = Sub(y, center_y);
        const auto dz      = Sub(z, center_z);
        const auto dist    = MulAdd(dx, dx, MulAdd(dy, dy, Mul(dz, dz)));
        const auto penalty = Mul(dist, penalty_vec);
        elect_total_trilinear += ReduceSum(d, IfThenElseZero(outside, penalty));
        emap_total_trilinear += ReduceSum(d, IfThenElseZero(outside, penalty));
      }

      // For atoms inside
      const auto inside = And(Not(outside), valid_coord);

      if (FindFirstTrue(d, inside) != -1) {
        // Calculate trilinear interpolation coordinates for in-bounds atoms
        x = Mul(Sub(x, min_x), inv_spacing_vec);
        y = Mul(Sub(y, min_y), inv_spacing_vec);
        z = Mul(Sub(z, min_z), inv_spacing_vec);

        // Map loading
        const auto atom_map_offsets = LoadN(di, map_ligand_offsets + index, remaining);

        // Decompose coordinates and weights
        const auto u0 = ConvertTo(di, x);
        const auto v0 = ConvertTo(di, y);
        const auto w0 = ConvertTo(di, z);

        const auto p0u = Sub(x, ConvertTo(d, u0));
        const auto p0v = Sub(y, ConvertTo(d, v0));
        const auto p0w = Sub(z, ConvertTo(d, w0));

        const auto p1u = Sub(one, p0u);
        const auto p1v = Sub(one, p0v);
        const auto p1w = Sub(one, p0w);

        // Precompute flattened indices
        const auto base_default_index =
            Add(Add(Add(Mul(v0, map_index_x_vec), Mul(w0, map_index_xy_vec)), u0), default_offset);
        elect_total_trilinear += ReduceSum(
            d,
            IfThenElseZero(inside,
                           Mul(trilinear_interpolation_vectorized<decltype(inside),
                                                                  decltype(p1u),
                                                                  decltype(base_default_index),
                                                                  fp_type>(electro_map,
                                                                           base_default_index,
                                                                           inside,
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
                               atom_charge)));
        dmap_total_trilinear += ReduceSum(
            d,
            IfThenElseZero(inside,
                           Mul(trilinear_interpolation_vectorized<decltype(inside),
                                                                  decltype(p1u),
                                                                  decltype(base_default_index),
                                                                  fp_type>(desolv_map,
                                                                           base_default_index,
                                                                           inside,
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
                               Abs(atom_charge))));
        const auto base_atom_index =
            Add(Add(Add(Mul(v0, map_index_x_vec), Mul(w0, map_index_xy_vec)), u0), atom_map_offsets);
        emap_total_trilinear += ReduceSum(
            d,
            IfThenElseZero(inside,
                           trilinear_interpolation_vectorized<decltype(inside),
                                                              decltype(p1u),
                                                              decltype(base_atom_index),
                                                              fp_type>(grid_maps,
                                                                       base_atom_index,
                                                                       inside,
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
                                                                       map_index_x_xy_plus_one_vec)));
      }
    }

    fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
    if (n_torsions > 0) {
      // TODO check how they are managed in memory these constants
      const auto ms_A_vec        = Set(d, mehler_solmajer::A);
      const auto ms_B_vec        = Set(d, mehler_solmajer::B);
      const auto ms_rk_vec       = Set(d, mehler_solmajer::rk);
      const auto ms_lambda_B_vec = Set(d, mehler_solmajer::lambda_B);
      // const auto ms_min_epsilon_vec = Set(d, std::numeric_limits<fp_type>::epsilon());
      const auto one_vec         = Set(d, fp_type{1});
      const auto elec_scale      = Set(d, ELECSCALE);
      const auto coeff_estat_vec = Set(d, autodock_parameters::coeff_estat);
      const auto nbc2_vec        = Set(d, nbc2);
      // const auto one_fp_vec         = Set(d, 1);
      const auto half_fp_vec = Set(d, -0.5);

      const auto xA_vec                      = Set(di, 12);
      const auto eintclamp_vec               = Set(d, EINTCLAMP);
      const auto rmin_elec_square_vec        = Set(d, RMIN_ELEC_SQUARE);
      const auto qsolpar_vec                 = Set(d, 0.01097);
      const auto coeff_desolv_vec            = Set(d, autodock_parameters::coeff_desolv);
      const auto reciprocal_sigma_square_vec = ApproximateReciprocal(Set(d, sigma_square));

      for (size_t i = 0; i < num_non_bond_loops * Lanes(d); i += Lanes(d)) {
        const auto remaining     = num_nonbond - i;
        const auto valid_nonbond = FirstN(d, remaining);

        const auto a1 = LoadN(di, non_bond_list_a1 + i, remaining);
        const auto a2 = LoadN(di, non_bond_list_a2 + i, remaining);

        // Load vectorized data for each coordinate difference between a1 and a2
        const auto diff_x =
            Sub(GatherIndexN(d, ligand_x, a1, remaining), GatherIndexN(d, ligand_x, a2, remaining));
        const auto diff_y =
            Sub(GatherIndexN(d, ligand_y, a1, remaining), GatherIndexN(d, ligand_y, a2, remaining));
        const auto diff_z =
            Sub(GatherIndexN(d, ligand_z, a1, remaining), GatherIndexN(d, ligand_z, a2, remaining));

        // Calculate squared distances
        const auto distance_two = Add(Mul(diff_x, diff_x), Add(Mul(diff_y, diff_y), Mul(diff_z, diff_z)));
        // Clamp `distance_two` between `RMIN_ELEC_SQUARE` and `distance_two` itself
        const auto clamped_distance_two = Max(rmin_elec_square_vec, distance_two);
        // Compute the square root of the clamped distance
        const auto distance = Mul(clamped_distance_two, ApproximateReciprocalSqrt(clamped_distance_two));

        //  Calculate  Electrostatic  Energy
        // TODO missing epsilon check
        // Calculate the r_dielectric term
        const auto epsilon = Add(
            Mul(ms_B_vec,
                ApproximateReciprocal(Add(one_vec, Mul(ms_rk_vec, Exp(d, Mul(ms_lambda_B_vec, distance)))))),
            ms_A_vec);
        // const auto ms_vec       = IfThenElse(Lt(epsilon, ms_min_epsilon_vec), one_fp_vec, epsilon);
        const auto r_dielectric = Mul(one_vec, ApproximateReciprocal(Mul(distance, epsilon)));

        // Calculate the electrostatic energy (e_elec)
        const auto charge_a1 = GatherIndexN(d, ligand_charge, a1, remaining);
        const auto charge_a2 = GatherIndexN(d, ligand_charge, a2, remaining);

        const auto e_elec =
            Mul(Mul(charge_a1, charge_a2), Mul(elec_scale, Mul(coeff_estat_vec, r_dielectric)));

        // Accumulate the electrostatic energy
        elect_total_eintcal += ReduceSum(d, IfThenElseZero(valid_nonbond, e_elec));

        // Calcuate desolv
        const auto ligand_vol_a1_v    = GatherIndexN(d, ligand_vol, a1, remaining);
        const auto ligand_vol_a2_v    = GatherIndexN(d, ligand_vol, a2, remaining);
        const auto ligand_solpar_a1_v = GatherIndexN(d, ligand_solpar, a1, remaining);
        const auto ligand_solpar_a2_v = GatherIndexN(d, ligand_solpar, a2, remaining);

        const auto nb_desolv =
            MulAdd(ligand_vol_a2_v,
                   Add(ligand_solpar_a1_v, Mul(qsolpar_vec, Abs(charge_a1))),
                   Mul(ligand_vol_a1_v, Add(ligand_solpar_a2_v, Mul(qsolpar_vec, Abs(charge_a2)))));

        const auto e_desolv = Mul(
            coeff_desolv_vec,
            Mul(Exp(d, Mul(Mul(half_fp_vec, reciprocal_sigma_square_vec), clamped_distance_two)), nb_desolv));
        dmap_total_eintcal += ReduceSum(d, IfThenElseZero(valid_nonbond, e_desolv));

        const auto low_distance = And(Lt(clamped_distance_two, nbc2_vec), valid_nonbond);
        auto e_vdW_Hb           = Zero(d);
        if (FindFirstTrue(d, low_distance) != -1) {
          const auto xB_vec = LoadN(di, xB_list + i, remaining);

          const auto xab_cond = Ne(xA_vec, xB_vec);

          if (FindFirstTrue(di, xab_cond) != -1) {
            const auto cA = LoadN(d, cA_list + i, remaining);
            const auto cB = LoadN(d, cB_list + i, remaining);

            const auto log_distance = Log(d, distance);

            const auto rA = Exp(d, Mul(ConvertTo(d, xA_vec), log_distance));
            const auto rB = Exp(d, Mul(ConvertTo(d, xB_vec), log_distance));

            // TODO IfThenElse to enable final reduction
            e_vdW_Hb = IfThenElseZero(
                RebindMask(d, xab_cond),
                Min(eintclamp_vec,
                    MulSub(cA, ApproximateReciprocal(rA), Mul(cB, ApproximateReciprocal(rB)))));
          }
        }
        emap_total_eintcal += ReduceSum(d, IfThenElseZero(low_distance, e_vdW_Hb));
      }
    }
    const fp_type tors_free_energy = n_torsions * autodock_parameters::coeff_tors;

    const fp_type total_trilinear = emap_total_trilinear + elect_total_trilinear + dmap_total_trilinear;
    const fp_type total_eintcal   = emap_total_eintcal + elect_total_eintcal + dmap_total_eintcal;
    return total_trilinear + total_eintcal + tors_free_energy;
  }

} // namespace mudock
