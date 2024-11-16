#include "mudock/type_alias.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <hwy/contrib/math/math-inl.h>
#include <hwy/highway.h>
#include <memory>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/ligand_maps.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/fjapp_utils.hpp>
#include <mudock/g_highway_implementation/geometric_transformations.hpp>
#include <mudock/g_highway_implementation/mutate.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/likwid_utils.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/utils.hpp>
#include <random>

#define FLATTENED_2D(x, y, index_x) ((y) *index_x + (x))

namespace mudock {
  // using HWY_NAMESPACE::translate_molecule;
  using namespace hwy::HWY_NAMESPACE;

  static constexpr auto coordinate_step = fp_type{0.2};
  static constexpr auto angle_step      = fp_type{4};

  // // TODO fix me from reorder_buffer.hpp
  // static constexpr int get_num_atom_clusters() { return 7; };
  // // the description of how we generate the clusters
  // static constexpr std::array<int, 7> atoms_clusters = {{0, 32, 64, 128, 160, 192, 256}};

  template<typename T>
  [[nodiscard]] inline const T random_gen_cpp(std::mt19937& generator,
                                              std::uniform_real_distribution<fp_type>& dist,
                                              const T& min,
                                              const T& max) {
    fp_type value;
    if constexpr (is_debug())
      value = fp_type{0.4};
    else {
      value = dist(generator);
    }
    return static_cast<T>(value * (max - min) + min);
  }

  inline int get_selection_distribution(std::mt19937& generator,
                                        std::uniform_real_distribution<fp_type>& dist,
                                        const int& population_number) {
    return random_gen_cpp<int>(generator, dist, 0, population_number - 1);
  };
  inline fp_type get_init_change_distribution(std::mt19937& generator,
                                              std::uniform_real_distribution<fp_type>& dist) {
    return random_gen_cpp<fp_type>(generator, dist, -45, 45);
  }
  inline fp_type get_mutation_change_distribution(std::mt19937& generator,
                                                  std::uniform_real_distribution<fp_type>& dist) {
    return random_gen_cpp<fp_type>(generator, dist, -10, 10);
  };
  inline fp_type get_mutation_coin_distribution(std::mt19937& generator,
                                                std::uniform_real_distribution<fp_type>& dist) {
    return random_gen_cpp<fp_type>(generator, dist, 0, 1);
  };
  inline int get_crossover_distribution(std::mt19937& generator,
                                        std::uniform_real_distribution<fp_type>& dist,
                                        const int& num_rotamers) {
    return random_gen_cpp<int>(generator, dist, 0, 6 + num_rotamers);
  };
  inline const chromosome& tournament_selection(std::mt19937& generator,
                                                std::uniform_real_distribution<fp_type>& dist,
                                                const int& tournament_length,
                                                const individual* __restrict__ population,
                                                const int& population_number) {
    const auto num_iterations = tournament_length;
    auto best_individual      = get_selection_distribution(generator, dist, population_number);
    for (int i = 0; i < num_iterations; ++i) {
      auto contendent = get_selection_distribution(generator, dist, population_number);
      if (population[contendent].score < population[best_individual].score) {
        best_individual = contendent;
      }
    }
    return population[best_individual].genes;
  }

  template<typename V, typename VI, typename T>
  inline V trilinear_interpolation_vectorized(const T* __restrict__ map,
                                              const VI& base_index,
                                              const V& p1u,
                                              const V& p1v,
                                              const V& p1w,
                                              const V& p0u,
                                              const V& p0v,
                                              const V& p0w,
                                              const int& map_index_x,
                                              const int& map_index_xy) {
    const HWY_FULL(T) d;
    const HWY_FULL(int) di;

    // Precompute flattened indices
    auto value = Zero(d);

    // Gather and accumulate the values based on offsets
    value = MulAdd(p1u * p1v * p1w, GatherIndex(d, map, base_index), value);
    value = MulAdd(p1u * p1v * p0w, GatherIndex(d, map, Add(base_index, Set(di, map_index_xy))), value);
    value = MulAdd(p1u * p0v * p1w, GatherIndex(d, map, Add(base_index, Set(di, map_index_x))), value);
    value = MulAdd(p1u * p0v * p0w,
                   GatherIndex(d, map, Add(base_index, Set(di, map_index_x + map_index_xy))),
                   value);
    value = MulAdd(p0u * p1v * p1w, GatherIndex(d, map, Add(base_index, Set(di, 1))), value);
    value = MulAdd(p0u * p1v * p0w, GatherIndex(d, map, Add(base_index, Set(di, 1 + map_index_xy))), value);
    value = MulAdd(p0u * p0v * p1w, GatherIndex(d, map, Add(base_index, Set(di, 1 + map_index_x))), value);
    value = MulAdd(p0u * p0v * p0w,
                   GatherIndex(d, map, Add(base_index, Set(di, 1 + map_index_x + map_index_xy))),
                   value);

    return value;
  }

  inline fp_type calc_energy(const fp_type* __restrict__ ligand_x,
                             const fp_type* __restrict__ ligand_y,
                             const fp_type* __restrict__ ligand_z,
                             const fp_type* __restrict__ ligand_vol,
                             const fp_type* __restrict__ ligand_solpar,
                             const fp_type* __restrict__ ligand_charge,
                             const int* __restrict__ ligand_num_hbond,
                             const fp_type* __restrict__ ligand_Rij_hb,
                             const fp_type* __restrict__ ligand_Rii,
                             const fp_type* __restrict__ ligand_epsij_hb,
                             const fp_type* __restrict__ ligand_epsii,
                             const int* __restrict__ map_ligand_offsets,
                             const int num_atoms,
                             const int n_torsions,
                             const int num_nonbond,
                             const int* __restrict__ non_bond_list_a1,
                             const int* __restrict__ non_bond_list_a2,
                             const fp_type* __restrict__ minimum,
                             const fp_type* __restrict__ maximum,
                             const fp_type* __restrict__ center,
                             const int map_index_x,
                             const int map_index_xy,
                             const fp_type* __restrict__ grid_maps,
                             const fp_type* __restrict__ electro_map,
                             const fp_type* __restrict__ desolv_map) {
    const HWY_FULL(fp_type) d;
    const HWY_FULL(int) di;

    fp_type elect_total_trilinear = 0;
    fp_type emap_total_trilinear  = 0;
    fp_type dmap_total_trilinear  = 0;

    for (int index = 0; index < num_atoms; index += Lanes(d)) {
      const auto remaining = num_atoms - index;
      // Load the x, y, z coordinates in a SIMD fashion
      auto x                 = LoadN(d, ligand_x + index, remaining);
      auto y                 = LoadN(d, ligand_y + index, remaining);
      auto z                 = LoadN(d, ligand_z + index, remaining);
      const auto valid_coord = FirstN(d, remaining);
      const auto atom_charge = LoadN(d, ligand_charge + index, remaining);

      // Set up constants for SIMD operations
      const auto min_x = Set(d, minimum[0]);
      const auto max_x = Set(d, maximum[0]);
      const auto min_y = Set(d, minimum[1]);
      const auto max_y = Set(d, maximum[1]);
      const auto min_z = Set(d, minimum[2]);
      const auto max_z = Set(d, maximum[2]);

      // Bounds check for each atom in vectorized form
      const auto outside_x = Or(Lt(x, min_x), Gt(x, max_x));
      const auto outside_y = Or(Lt(y, min_y), Gt(y, max_y));
      const auto outside_z = Or(Lt(z, min_z), Gt(z, max_z));
      const auto outside   = And(Or(outside_x, Or(outside_y, outside_z)), valid_coord);

      // For atoms outside the boundaries
      if (FindFirstTrue(d, outside) != -1) {
        const auto dx      = Sub(x, Set(d, center[0]));
        const auto dy      = Sub(y, Set(d, center[1]));
        const auto dz      = Sub(z, Set(d, center[2]));
        const auto dist    = MulAdd(dx, dx, MulAdd(dy, dy, Mul(dz, dz)));
        const auto penalty = Mul(dist, Set(d, ENERGYPENALTY));
        elect_total_trilinear += ReduceSum(d, IfThenElseZero(outside, penalty));
        emap_total_trilinear += ReduceSum(d, IfThenElseZero(outside, penalty));
      }

      // For atoms inside
      const auto inside = And(Not(outside), valid_coord);

      if (FindFirstTrue(d, inside) != -1) {
        // Calculate trilinear interpolation coordinates for in-bounds atoms
        x = Mul(Sub(x, min_x), Set(d, inv_spacing));
        y = Mul(Sub(y, min_y), Set(d, inv_spacing));
        z = Mul(Sub(z, min_z), Set(d, inv_spacing));

        // Map loading
        const auto atom_map_offsets = LoadN(di, map_ligand_offsets + index, remaining);

        // Decompose coordinates and weights
        const auto u0 = ConvertTo(di, x);
        const auto v0 = ConvertTo(di, y);
        const auto w0 = ConvertTo(di, z);

        const auto p0u = Sub(x, ConvertTo(d, u0));
        const auto p0v = Sub(y, ConvertTo(d, v0));
        const auto p0w = Sub(z, ConvertTo(d, w0));

        const auto one = Set(d, fp_type{1.0});
        const auto p1u = Sub(one, p0u);
        const auto p1v = Sub(one, p0v);
        const auto p1w = Sub(one, p0w);

        // Precompute flattened indices
        const auto default_offset = Set(di, 0);
        const auto base_default_index =
            Add(Add(Add(Mul(v0, Set(di, map_index_x)), Mul(w0, Set(di, map_index_xy))), u0), default_offset);
        elect_total_trilinear += ReduceSum(
            d,
            IfThenElseZero(
                inside,
                trilinear_interpolation_vectorized<decltype(p1u), decltype(base_default_index), fp_type>(
                    electro_map,
                    base_default_index,
                    p1u,
                    p1v,
                    p1w,
                    p0u,
                    p0v,
                    p0w,
                    map_index_x,
                    map_index_xy)) *
                atom_charge);
        dmap_total_trilinear += ReduceSum(
            d,
            IfThenElseZero(
                inside,
                trilinear_interpolation_vectorized<decltype(p1u), decltype(base_default_index), fp_type>(
                    desolv_map,
                    base_default_index,
                    p1u,
                    p1v,
                    p1w,
                    p0u,
                    p0v,
                    p0w,
                    map_index_x,
                    map_index_xy)) *
                Abs(atom_charge));
        const auto base_atom_index =
            Add(Add(Add(Mul(v0, Set(di, map_index_x)), Mul(w0, Set(di, map_index_xy))), u0),
                atom_map_offsets);
        emap_total_trilinear += ReduceSum(
            d,
            IfThenElseZero(
                inside,
                trilinear_interpolation_vectorized<decltype(p1u), decltype(base_atom_index), fp_type>(
                    grid_maps,
                    base_atom_index,
                    p1u,
                    p1v,
                    p1w,
                    p0u,
                    p0v,
                    p0w,
                    map_index_x,
                    map_index_xy)));
      }
    }

    fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
    if (n_torsions > 0) {
      const auto ms_A_vec           = Set(d, mehler_solmajer::A);
      const auto ms_B_vec           = Set(d, mehler_solmajer::B);
      const auto ms_rk_vec          = Set(d, mehler_solmajer::rk);
      const auto ms_lambda_B_vec    = Set(d, mehler_solmajer::lambda_B);
      const auto ms_min_epsilon_vec = Set(d, std::numeric_limits<fp_type>::epsilon());
      const auto one_vec            = Set(d, fp_type{1});
      const auto elec_scale         = Set(d, ELECSCALE);
      const auto coeff_estat_vec    = Set(d, autodock_parameters::coeff_estat);
      const auto nbc2_vec           = Set(d, nbc2);
      const auto hbond_two_vev      = Set(di, 2);
      const auto two_fp_vec         = Set(d, 2);
      const auto one_fp_vec         = Set(d, 1);
      const auto half_fp_vec        = Set(d, -0.5);

      const auto hbond_one_vec    = Set(di, 1);
      const auto xA_vec           = Set(di, 12);
      const auto xB_ten_vec       = Set(di, 10);
      const auto xB_six_vec       = Set(di, 6);
      const auto eintclamp_vec    = Set(d, EINTCLAMP);
      const auto qsolpar_vec      = Set(d, 0.01097);
      const auto coeff_desolv_vec = Set(d, autodock_parameters::coeff_desolv);
      const auto sigma_square_vec = Set(d, sigma_square);

      for (int i = 0; i < num_nonbond; i += Lanes(di)) {
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
        const auto clamped_distance_two = Max(Set(d, RMIN_ELEC_SQUARE), distance_two);
        // Compute the square root of the clamped distance
        const auto distance = Sqrt(clamped_distance_two);

        //  Calculate  Electrostatic  Energy
        // TODO missing epsilon check
        // Calculate the r_dielectric term
        const auto epsilon =
            Add(Div(ms_B_vec, Add(one_vec, Mul(ms_rk_vec, Exp(d, Mul(ms_lambda_B_vec, distance))))),
                ms_A_vec);
        const auto ms_vec       = IfThenElse(Lt(epsilon, ms_min_epsilon_vec), one_fp_vec, epsilon);
        const auto r_dielectric = Div(one_vec, Mul(distance, ms_vec));

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
            Add(Mul(ligand_vol_a2_v, Add(ligand_solpar_a1_v, Mul(qsolpar_vec, Abs(charge_a1)))),
                Mul(ligand_vol_a1_v, Add(ligand_solpar_a2_v, Mul(qsolpar_vec, Abs(charge_a2)))));

        const auto e_desolv =
            Mul(coeff_desolv_vec,
                Mul(Exp(d, Mul(Div(half_fp_vec, sigma_square_vec), clamped_distance_two)), nb_desolv));
        dmap_total_eintcal += ReduceSum(d, IfThenElseZero(valid_nonbond, e_desolv));

        const auto low_distance = And(Lt(clamped_distance_two, nbc2_vec), valid_nonbond);
        auto e_vdW_Hb           = Zero(d);
        if (FindFirstTrue(d, low_distance) != -1) {
          // Vectorized arrays
          const auto hbond_i_v = GatherIndexN(di, ligand_num_hbond, a1, remaining);
          const auto hbond_j_v =
              GatherIndexN(di, ligand_num_hbond, a2, remaining); // example: next element for hbond_j

          const auto Rij_hb_i_v = GatherIndexN(d, ligand_Rij_hb, a1, remaining);
          const auto Rij_hb_j_v = GatherIndexN(d, ligand_Rij_hb, a2, remaining);

          const auto Rii_i_v = GatherIndexN(d, ligand_Rii, a1, remaining);
          const auto Rii_j_v = GatherIndexN(d, ligand_Rii, a2, remaining);

          const auto epsij_hb_i_v = GatherIndexN(d, ligand_epsij_hb, a1, remaining);
          const auto epsij_hb_j_v = GatherIndexN(d, ligand_epsij_hb, a2, remaining);

          const auto epsii_i_v = GatherIndexN(d, ligand_epsii, a1, remaining);
          const auto epsii_j_v = GatherIndexN(d, ligand_epsii, a2, remaining);

          // Compute Rij and epsij with masks
          // TODO force to do floating point comparison due to IfThenElse
          const auto i_donor_j_acceptor =
              RebindMask(d,
                         And(Or(Eq(hbond_i_v, hbond_one_vec), Eq(hbond_i_v, hbond_two_vev)),
                             Gt(hbond_j_v, hbond_two_vev)));
          const auto j_donor_i_acceptor =
              RebindMask(d,
                         And(Or(Eq(hbond_j_v, hbond_one_vec), Eq(hbond_j_v, hbond_two_vev)),
                             Gt(hbond_i_v, hbond_two_vev)));

          const auto Rij_v =
              IfThenElse(i_donor_j_acceptor,
                         Rij_hb_j_v,
                         IfThenElse(j_donor_i_acceptor, Rij_hb_i_v, (Rii_i_v + Rii_j_v) / two_fp_vec));

          const auto epsij_v =
              IfThenElse(i_donor_j_acceptor,
                         epsij_hb_j_v,
                         IfThenElse(j_donor_i_acceptor, epsij_hb_i_v, Sqrt(Mul(epsii_i_v, epsii_j_v))));

          const auto xB_vec =
              IfThenElse(RebindMask(di, Or(i_donor_j_acceptor, j_donor_i_acceptor)), xB_ten_vec, xB_six_vec);

          const auto xab_cond = Ne(xA_vec, xB_vec);

          if (FindFirstTrue(di, xab_cond) != -1) {
            const auto tmp = Div(epsij_v, ConvertTo(d, Sub(xA_vec, xB_vec)));

            // Hack, Pow is not yet support by GH [x^y = exp(y*ln(x))]
            const auto pow_A_vec = Exp(d, Mul(ConvertTo(d, xA_vec), Log(d, Rij_v)));
            const auto pow_B_vec = Exp(d, Mul(ConvertTo(d, xB_vec), Log(d, Rij_v)));
            const auto cA        = Mul(tmp, Mul(pow_A_vec, ConvertTo(d, xB_vec)));
            const auto cB        = Mul(tmp, Mul(pow_B_vec, ConvertTo(d, xA_vec)));

            const auto rA = Exp(d, Mul(ConvertTo(d, xA_vec), Log(d, distance)));
            const auto rB = Exp(d, Mul(ConvertTo(d, xB_vec), Log(d, distance)));

            // TODO IfThenElse to enable final reduction
            e_vdW_Hb =
                IfThenElseZero(RebindMask(d, xab_cond), Min(eintclamp_vec, Sub(Div(cA, rA), Div(cB, rB))));
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

  void evaluate_fitness(const fp_type* __restrict__ ligand_x,
                        const fp_type* __restrict__ ligand_y,
                        const fp_type* __restrict__ ligand_z,
                        const fp_type* __restrict__ ligand_vol,
                        const fp_type* __restrict__ ligand_solpar,
                        const fp_type* __restrict__ ligand_charge,
                        const int* __restrict__ ligand_num_hbond,
                        const fp_type* __restrict__ ligand_Rij_hb,
                        const fp_type* __restrict__ ligand_Rii,
                        const fp_type* __restrict__ ligand_epsij_hb,
                        const fp_type* __restrict__ ligand_epsii,
                        const int* __restrict__ map_ligand_offsets,
                        const int num_atoms,
                        const int num_rotamers,
                        const int* __restrict__ frag_masks,
                        const int* __restrict__ frag_start_indexes,
                        const int* __restrict__ frag_stop_indexes,
                        const int num_nonbond,
                        const int* __restrict__ non_bond_list_a1,
                        const int* __restrict__ non_bond_list_a2,
                        const fp_type* __restrict__ grid_maps,
                        const fp_type* __restrict__ electro_map,
                        const fp_type* __restrict__ desolv_map,
                        const int num_generations,
                        const int population_size,
                        const int tournament_length,
                        const fp_type mutation_prob,
                        const fp_type* __restrict__ minimum,
                        const fp_type* __restrict__ maximum,
                        const fp_type* __restrict__ center,
                        const int map_index_x,
                        const int map_index_xy,
                        individual* __restrict__ population_buffer1,
                        individual* __restrict__ population_buffer2,
                        const int seed) {
    std::uniform_real_distribution<fp_type> dist{fp_type{0.0}, fp_type{1.0}};
    std::mt19937 generator(seed);

    auto* population      = population_buffer1;
    auto* next_population = population_buffer2;
    auto altered_x        = std::make_unique<std::array<fp_type, max_static_atoms()>>();
    auto altered_y        = std::make_unique<std::array<fp_type, max_static_atoms()>>();
    auto altered_z        = std::make_unique<std::array<fp_type, max_static_atoms()>>();

    FJAPP_MARKER_START("GA");
    LIKWID_MARKER_START("GA");

// Randomly initialize the population
#pragma clang loop unroll(enable)
    for (int element_index = 0; element_index < population_size; ++element_index) {
      auto& element = population[element_index];
#pragma clang loop unroll(enable)
      for (int i{0}; i < 3; ++i) { // initialize the rigid translation
        element.genes[i] = get_init_change_distribution(generator, dist) * coordinate_step;
      }
#pragma clang loop unroll(enable)
      for (int i{3}; i < 6 + num_rotamers; ++i) { // initialize the rotations
        element.genes[i] = get_init_change_distribution(generator, dist) * angle_step;
      }
    }

    for (int generation = 0; generation < num_generations; ++generation) {
      // Evaluate the fitness of the population
      for (int element_index = 0; element_index < population_size; ++element_index) {
        auto& element = population[element_index];

        std::memcpy(altered_x.get()->data(), ligand_x, num_atoms * sizeof(fp_type));
        std::memcpy(altered_y.get()->data(), ligand_y, num_atoms * sizeof(fp_type));
        std::memcpy(altered_z.get()->data(), ligand_z, num_atoms * sizeof(fp_type));

        // TODO check it it makes sense -> print the MOL2
        // apply the transformation encoded in the element genes to the original ligand
        apply(altered_x.get()->data(),
              altered_y.get()->data(),
              altered_z.get()->data(),
              element.genes,
              num_atoms,
              num_rotamers,
              frag_masks,
              frag_start_indexes,
              frag_stop_indexes);

        // compute the energy of the system
        const auto energy = calc_energy(altered_x.get()->data(),
                                        altered_y.get()->data(),
                                        altered_z.get()->data(),
                                        ligand_vol,
                                        ligand_solpar,
                                        ligand_charge,
                                        ligand_num_hbond,
                                        ligand_Rij_hb,
                                        ligand_Rii,
                                        ligand_epsij_hb,
                                        ligand_epsii,
                                        map_ligand_offsets,
                                        num_atoms,
                                        num_rotamers,
                                        num_nonbond,
                                        non_bond_list_a1,
                                        non_bond_list_a2,
                                        minimum,
                                        maximum,
                                        center,
                                        map_index_x,
                                        map_index_xy,
                                        grid_maps,
                                        electro_map,
                                        desolv_map);
        element.score     = energy; // dummy implementation to test the genetic
      }

      // Generate the new population
      for (int element_index = 0; element_index < population_size; ++element_index) {
        auto& next_individual = next_population[element_index];
        // select the parent
        const auto& parent1 =
            tournament_selection(generator, dist, tournament_length, population, population_size);
        const auto& parent2 =
            tournament_selection(generator, dist, tournament_length, population, population_size);

        // generate the offspring
        const auto split_index = get_crossover_distribution(generator, dist, num_rotamers);
        std::copy(std::begin(parent1), std::begin(parent1) + split_index, std::begin(next_individual.genes));
        std::copy(std::begin(parent2) + split_index,
                  std::end(parent2),
                  std::begin(next_individual.genes) + split_index);
        next_individual.score = fp_type{0};

// mutate the offspring
#pragma clang loop unroll(enable)
        for (int i{0}; i < 3; ++i) {
          if (get_mutation_coin_distribution(generator, dist) < mutation_prob)
            next_individual.genes[i] += get_mutation_change_distribution(generator, dist) * coordinate_step;
        }
#pragma clang loop unroll(enable)
        for (int i{3}; i < 6 + num_rotamers; ++i) {
          if (get_mutation_coin_distribution(generator, dist) < mutation_prob)
            next_individual.genes[i] += get_mutation_change_distribution(generator, dist) * angle_step;
        }
      }

      // swap the new population with the old one
      const auto temp = population;
      population      = next_population;
      next_population = temp;
    }
    FJAPP_MARKER_STOP("GA");
    LIKWID_MARKER_STOP("GA");
  }
} // namespace mudock
