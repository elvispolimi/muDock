#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/ligand_maps.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/fjapp_utils.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/likwid_utils.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/utils.hpp>
#include <random>

#define FLATTENED_2D(x, y, index_x)              ((y) * index_x + (x))
#define FLATTENED_3D(x, y, z, index_x, index_xy) (index_xy * (z) + (y) * index_x + (x))

namespace mudock {
  static constexpr auto coordinate_step = fp_type{0.2};
  static constexpr auto angle_step      = fp_type{4};

  // TODO fix me from reorder_buffer.hpp
  static constexpr int get_num_atom_clusters() { return 7; };
  // the description of how we generate the clusters
  static constexpr std::array<int, 7> atoms_clusters = {{0, 32, 64, 128, 160, 192, 256}};

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

  inline fp_type trilinear_interpolation(const fp_type* __restrict__ map,
                                         const fp_type* __restrict__ coeffs,
                                         const int base_index,
                                         const int& map_index_x,
                                         const int& map_index_xy) {
    fp_type value{0};

    value = coeffs[0] * map[base_index] + value;
    value = coeffs[1] * map[base_index + map_index_xy] + value;
    value = coeffs[2] * map[base_index + map_index_x] + value;
    value = coeffs[3] * map[base_index + map_index_x + map_index_xy] + value;
    value = coeffs[4] * map[base_index + 1] + value;
    value = coeffs[5] * map[base_index + 1 + map_index_xy] + value;
    value = coeffs[6] * map[base_index + 1 + map_index_x] + value;
    value = coeffs[7] * map[base_index + 1 + map_index_x + map_index_xy] + value;

    return value;
  }

  template<int NUM_ATOMS>
  inline void translate_molecule(fp_type* __restrict__ x,
                                 fp_type* __restrict__ y,
                                 fp_type* __restrict__ z,
                                 const int num_atoms,
                                 const fp_type offset_x,
                                 const fp_type offset_y,
                                 const fp_type offset_z) {
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
#pragma clang loop unroll(enable)
    for (int i = 0; i < NUM_ATOMS; ++i)
      if (i < num_atoms) {
        x[i] += offset_x;
        y[i] += offset_y;
        z[i] += offset_z;
      }
  }

  template<int NUM_ATOMS>
  inline void rotate_molecule(fp_type* __restrict__ x,
                              fp_type* __restrict__ y,
                              fp_type* __restrict__ z,
                              const int num_atoms,
                              const fp_type angle_x,
                              const fp_type angle_y,
                              const fp_type angle_z) {
    // compute the molecule center of mass
    point3D c{0, 0, 0};
#pragma GCC ivdep
#pragma clang loop vectorize(enable)
#pragma clang loop unroll(enable)
    for (int i = 0; i < NUM_ATOMS; i++)
      if (i < num_atoms) {
        c.x += x[i];
        c.y += y[i];
        c.z += z[i];
      }
    c.x /= num_atoms;
    c.y /= num_atoms;
    c.z /= num_atoms;

    // compute the angles sine and cosine
    const auto rad_x = deg_to_rad(angle_x), rad_y = deg_to_rad(angle_y), rad_z = deg_to_rad(angle_z);
    const auto cx = std::cos(rad_x), sx = std::sin(rad_x);
    const auto cy = std::cos(rad_y), sy = std::sin(rad_y);
    const auto cz = std::cos(rad_z), sz = std::sin(rad_z);

    // compute the rotation matrix defined as Rz*Ry*Rx
    const auto m00 = cy * cz;
    const auto m01 = sx * sy * cz - cx * sz;
    const auto m02 = cx * sy * cz + sx * sz;
    const auto m10 = cy * sz;
    const auto m11 = sx * sy * sz + cx * cz;
    const auto m12 = cx * sy * sz - sx * cz;
    const auto m20 = -sy;
    const auto m21 = sx * cy;
    const auto m22 = cx * cy;

// apply the rotation matrix
#pragma GCC ivdep
#pragma clang loop vectorize(enable)
#pragma clang loop unroll(enable)
    for (int i = 0; i < NUM_ATOMS; ++i)
      if (i < num_atoms) {
        const auto translated_x = x[i] - c.x, translated_y = y[i] - c.y, translated_z = z[i] - c.z;
        x[i] = translated_x * m00 + translated_y * m01 + translated_z * m02 + c.x;
        y[i] = translated_x * m10 + translated_y * m11 + translated_z * m12 + c.y;
        z[i] = translated_x * m20 + translated_y * m21 + translated_z * m22 + c.z;
      }
  }

  template<int NUM_ATOMS>
  inline void rotate_fragment(fp_type* __restrict__ x,
                              fp_type* __restrict__ y,
                              fp_type* __restrict__ z,
                              const int num_atoms,
                              const int* __restrict__ frag_mask,
                              const int start_index,
                              const int stop_index,
                              const fp_type angle) {
    // compute the axis vector (and some properties)
    const auto origx = x[start_index], origy = y[start_index], origz = z[start_index];
    const auto destx = x[stop_index], desty = y[stop_index], destz = z[stop_index];
    const auto u  = destx - origx;
    const auto v  = desty - origy;
    const auto w  = destz - origz;
    const auto u2 = u * u, v2 = v * v, w2 = w * w;
    const auto l2 = u * u + v * v + w * w;
    const auto l  = std::sqrt(l2);

    // compute the angle sine and cosine
    const auto rad = deg_to_rad(angle);
    const auto s = std::sin(rad), c = std::cos(rad);
    const auto one_minus_c = fp_type{1} - c;
    const auto ls          = l * s;

    // Precompute common sub-expressions to reduce redundant calculations
    const auto inv_l2 = fp_type{1} / l2;
    const auto us_vc  = u * v * one_minus_c;
    const auto uw_vc  = u * w * one_minus_c;
    const auto vw_vc  = v * w * one_minus_c;

    // compute the rotation matrix (rodrigues' rotation formula)
    const auto m00 = (u2 + (v2 + w2) * c) * inv_l2;
    const auto m01 = (us_vc - w * l * s) * inv_l2;
    const auto m02 = (uw_vc + v * l * s) * inv_l2;
    const auto m03 =
        ((origx * (v2 + w2) - u * (origy * v + origz * w)) * one_minus_c + (origy * w - origz * v) * ls) *
        inv_l2;

    const auto m10 = (us_vc + w * ls) * inv_l2;
    const auto m11 = (v2 + (u2 + w2) * c) * inv_l2;
    const auto m12 = (vw_vc - u * ls) * inv_l2;
    const auto m13 =
        ((origy * (u2 + w2) - v * (origx * u + origz * w)) * one_minus_c + (origz * u - origx * w) * ls) *
        inv_l2;

    const auto m20 = (uw_vc - v * ls) * inv_l2;
    const auto m21 = (vw_vc + u * ls) * inv_l2;
    const auto m22 = (w2 + (u2 + v2) * c) * inv_l2;
    const auto m23 =
        ((origz * (u2 + v2) - w * (origx * u + origy * v)) * one_minus_c + (origx * v - origy * u) * ls) *
        inv_l2;

// apply the rotation matrix
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
#pragma clang loop unroll(enable)
    for (int i = 0; i < NUM_ATOMS; ++i) {
      if (frag_mask[i] != 0 && i < num_atoms) {
        const auto prev_x = x[i], prev_y = y[i], prev_z = z[i];
        x[i] = prev_x * m00 + prev_y * m01 + prev_z * m02 + m03;
        y[i] = prev_x * m10 + prev_y * m11 + prev_z * m12 + m13;
        z[i] = prev_x * m20 + prev_y * m21 + prev_z * m22 + m23;
      }
    }
  }

  template<int NUM_ATOMS>
  inline void apply(fp_type* __restrict__ x,
                    fp_type* __restrict__ y,
                    fp_type* __restrict__ z,
                    const chromosome& c,
                    const int num_atoms,
                    const int num_rotamers,
                    const int* __restrict__ frag_masks,
                    const int* __restrict__ frag_start_indexes,
                    const int* __restrict__ frag_stop_indexes) {
    // apply rigid transformations
    translate_molecule<NUM_ATOMS>(x, y, z, num_atoms, c[0], c[1], c[2]);
    rotate_molecule<NUM_ATOMS>(x, y, z, num_atoms, c[3], c[4], c[5]);

// change the molecule shape
#pragma clang loop unroll(enable)
    for (int i = 0; i < num_rotamers; ++i) {
      const auto* bitmask    = frag_masks + i * num_atoms;
      const auto start_index = frag_start_indexes[i];
      const auto stop_index  = frag_stop_indexes[i];
      rotate_fragment<NUM_ATOMS>(x, y, z, num_atoms, bitmask, start_index, stop_index, c[int{6} + i]);
    }
  }

  template<int NUM_ATOMS>
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
                             const ligand_map_types* __restrict__ map_ligand_types,
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
                             const fp_type* const __restrict__* const __restrict__ grid_maps,
                             const fp_type* __restrict__ electro_map,
                             const fp_type* __restrict__ desolv_map) {
    fp_type elect_total_trilinear = 0;
    fp_type emap_total_trilinear  = 0;
    fp_type dmap_total_trilinear  = 0;

#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
#pragma clang loop unroll(enable)
    for (int index = 0; index < NUM_ATOMS; ++index)
      if (index < num_atoms) {
        fp_type coord[3]{ligand_x[index], ligand_y[index], ligand_z[index]};

        if (coord[0] < minimum[0] || coord[0] > maximum[0] || coord[1] < minimum[1] ||
            coord[1] > maximum[1] || coord[2] < minimum[2] || coord[2] > maximum[2]) {
          const fp_type dist = std::pow(coord[0] - center[0], fp_type{2}) +
                               std::pow(coord[1] - center[1], fp_type{2}) +
                               std::pow(coord[2] - center[2], fp_type{2});
          const fp_type epenalty = dist * ENERGYPENALTY;
          elect_total_trilinear += epenalty;
          emap_total_trilinear += epenalty;
        } else {
          const auto& atom_charge = ligand_charge[index];
          const fp_type* atom_map = grid_maps[static_cast<int>(map_ligand_types[index])];

          coord[0] = (coord[0] - minimum[0]) * inv_spacing;
          coord[1] = (coord[1] - minimum[1]) * inv_spacing;
          coord[2] = (coord[2] - minimum[2]) * inv_spacing;

          const int u0      = coord[0];
          const fp_type p0u = coord[0] - static_cast<fp_type>(u0);
          const fp_type p1u = fp_type{1} - p0u;

          const int v0      = coord[1];
          const fp_type p0v = coord[1] - static_cast<fp_type>(v0);
          const fp_type p1v = fp_type{1} - p0v;

          const int w0      = coord[2];
          const fp_type p0w = coord[2] - static_cast<fp_type>(w0);
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

          // Precompute flattened indices
          const int base_index = FLATTENED_3D(u0, v0, w0, map_index_x, map_index_xy);
          // Trilinear Interpolationp
          elect_total_trilinear +=
              trilinear_interpolation(electro_map, coeffs, base_index, map_index_x, map_index_xy) *
              atom_charge;
          emap_total_trilinear +=
              trilinear_interpolation(atom_map, coeffs, base_index, map_index_x, map_index_xy);
          dmap_total_trilinear +=
              trilinear_interpolation(desolv_map, coeffs, base_index, map_index_x, map_index_xy) *
              std::fabs(atom_charge);
        }
      }

    fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
    if (n_torsions > 0) {
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
#pragma clang loop unroll(enable)
#pragma fj loop prefetch
#pragma statement scache_isolate_assign ligand_x, ligand_y, ligand_z, ligand_charge, ligand_num_hbond, \
    ligand_Rij_hb, ligand_Rii, ligand_epsij_hb, ligand_epsii
      for (int i = 0; i < num_nonbond; ++i) {
        const int& a1 = non_bond_list_a1[i];
        const int& a2 = non_bond_list_a2[i];

        const fp_type distance_two = std::pow(ligand_x[a1] - ligand_x[a2], fp_type{2}) +
                                     std::pow(ligand_y[a1] - ligand_y[a2], fp_type{2}) +
                                     std::pow(ligand_z[a1] - ligand_z[a2], fp_type{2});
        const fp_type distance_two_clamp = std::clamp(distance_two, RMIN_ELEC_SQUARE, distance_two);
        const fp_type distance           = std::sqrt(distance_two_clamp);

        //  Calculate  Electrostatic  Energy
        const fp_type r_dielectric = fp_type{1} / (distance * calc_ddd_Mehler_Solmajer(distance));
        const fp_type e_elec       = ligand_charge[a1] * ligand_charge[a2] * ELECSCALE *
                               autodock_parameters::coeff_estat * r_dielectric;
        elect_total_eintcal += e_elec;

        // Calcuare desolv
        const fp_type nb_desolv =
            (ligand_vol[a2] * (ligand_solpar[a1] + qsolpar * std::fabs(ligand_charge[a1])) +
             ligand_vol[a1] * (ligand_solpar[a2] + qsolpar * std::fabs(ligand_charge[a2])));

        const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                 std::exp(fp_type{-0.5} / (sigma_square) *distance_two_clamp) * nb_desolv;
        dmap_total_eintcal += e_desolv;
        fp_type e_vdW_Hb{0};
        if (distance_two_clamp < nbc2) {
          //  Find internal energy parameters, i.e.  epsilon and r-equilibrium values...
          //  Lennard-Jones and Hydrogen Bond Potentials
          // This can be precomputed as in intnbtable.cc
          const auto& hbond_i    = ligand_num_hbond[a1];
          const auto& hbond_j    = ligand_num_hbond[a2];
          const auto& Rij_hb_i   = ligand_Rij_hb[a1];
          const auto& Rij_hb_j   = ligand_Rij_hb[a2];
          const auto& Rii_i      = ligand_Rii[a1];
          const auto& Rii_j      = ligand_Rii[a2];
          const auto& epsij_hb_i = ligand_epsij_hb[a1];
          const auto& epsij_hb_j = ligand_epsij_hb[a2];
          const auto& epsii_i    = ligand_epsii[a1];
          const auto& epsii_j    = ligand_epsii[a2];

          // we need to determine the correct xA and xB exponents
          int xA = 12; // for both LJ, 12-6 and HB, 12-10, xA is 12
          int xB = 6;  // assume we have LJ, 12-6

          fp_type Rij{0}, epsij{0};
          if ((hbond_i == 1 || hbond_i == 2) && hbond_j > 2) {
            // i is a donor and j is an acceptor.
            // i is a hydrogen, j is a heteroatom
            Rij   = Rij_hb_j;
            epsij = epsij_hb_j;
            xB    = 10;
          } else if ((hbond_i > 2) && (hbond_j == 1 || hbond_j == 2)) {
            // i is an acceptor and j is a donor.
            // i is a heteroatom, j is a hydrogen
            Rij   = Rij_hb_i;
            epsij = epsij_hb_i;
            xB    = 10;
          } else {
            // we need to calculate the arithmetic mean of Ri and Rj
            Rij = (Rii_i + Rii_j) / fp_type{2};
            // we need to calculate the geometric mean of epsi and epsj
            epsij = std::sqrt(epsii_i * epsii_j);
          }
          if (xA != xB) {
            const fp_type tmp = epsij / (xA - xB);
            const fp_type cA  = tmp * std::pow(Rij, static_cast<fp_type>(xA)) * xB;
            const fp_type cB  = tmp * std::pow(Rij, static_cast<fp_type>(xB)) * xA;

            const fp_type rA = std::pow(distance, static_cast<fp_type>(xA));
            const fp_type rB = std::pow(distance, static_cast<fp_type>(xB));

            e_vdW_Hb = std::min(EINTCLAMP, (cA / rA - cB / rB));
          }
        }
        emap_total_eintcal += e_vdW_Hb;
      }
    }
    const fp_type tors_free_energy = n_torsions * autodock_parameters::coeff_tors;

    const fp_type total_trilinear = emap_total_trilinear + elect_total_trilinear + dmap_total_trilinear;
    const fp_type total_eintcal   = emap_total_eintcal + elect_total_eintcal + dmap_total_eintcal;
    return total_trilinear + total_eintcal + tors_free_energy;
  }

  template<int NUM_ATOMS>
  void evaluate_fitness_impl(const fp_type* __restrict__ ligand_x,
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
                             const ligand_map_types* __restrict__ map_ligand_types,
                             const int num_atoms,
                             const int num_rotamers,
                             const int* __restrict__ frag_masks,
                             const int* __restrict__ frag_start_indexes,
                             const int* __restrict__ frag_stop_indexes,
                             const int num_nonbond,
                             const int* __restrict__ non_bond_list_a1,
                             const int* __restrict__ non_bond_list_a2,
                             const fp_type* const __restrict__* const __restrict__ grid_maps,
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
    // TODO enable vectorization
    for (int element_index = 0; element_index < population_size; ++element_index) {
      auto& element = population[element_index];
#pragma clang loop unroll(enable)
      for (int i{0}; i < 3; ++i) { // initialize the rigid translation
        element.genes[i] = get_init_change_distribution(generator, dist) * coordinate_step;
      }
      for (int i{3}; i < 6 + num_rotamers; ++i) { // initialize the rotations
        element.genes[i] = get_init_change_distribution(generator, dist) * angle_step;
      }
    }

    for (int generation = 0; generation < num_generations; ++generation) {
// Evaluate the fitness of the population
#pragma fj loop prefetch
#pragma statement scache_isolate_assign ligand_x, ligand_y, ligand_z
      for (int element_index = 0; element_index < population_size; ++element_index) {
        auto& element = population[element_index];

        std::memcpy(altered_x.get()->data(), ligand_x, num_atoms * sizeof(fp_type));
        std::memcpy(altered_y.get()->data(), ligand_y, num_atoms * sizeof(fp_type));
        std::memcpy(altered_z.get()->data(), ligand_z, num_atoms * sizeof(fp_type));

        // TODO check it it makes sense -> print the MOL2
        // apply the transformation encoded in the element genes to the original ligand
        apply<NUM_ATOMS>(altered_x.get()->data(),
                         altered_y.get()->data(),
                         altered_z.get()->data(),
                         element.genes,
                         num_atoms,
                         num_rotamers,
                         frag_masks,
                         frag_start_indexes,
                         frag_stop_indexes);

        // compute the energy of the system
        const auto energy = calc_energy<NUM_ATOMS>(altered_x.get()->data(),
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
                                                   map_ligand_types,
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
      // TODO enable vectorization
      for (int element_index = 0; element_index < population_size; ++element_index) {
        auto& next_individual = next_population[element_index];
        // select the parent
        auto best_individual_1 = get_selection_distribution(generator, dist, population_size);
        auto best_individual_2 = get_selection_distribution(generator, dist, population_size);
        for (int i = 0; i < tournament_length; ++i) {
          const auto contendent_1 = get_selection_distribution(generator, dist, population_size);
          const auto contendent_2 = get_selection_distribution(generator, dist, population_size);
          if (population[contendent_1].score < population[best_individual_1].score) {
            best_individual_1 = contendent_1;
          }
          if (population[contendent_2].score < population[best_individual_2].score) {
            best_individual_2 = contendent_2;
          }
        }
        const auto& parent1 = population[best_individual_1].genes;
        const auto& parent2 = population[best_individual_2].genes;

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
                        const ligand_map_types* __restrict__ map_ligand_types,
                        const int num_atoms,
                        const int num_rotamers,
                        const int* __restrict__ frag_masks,
                        const int* __restrict__ frag_start_indexes,
                        const int* __restrict__ frag_stop_indexes,
                        const int num_nonbond,
                        const int* __restrict__ non_bond_list_a1,
                        const int* __restrict__ non_bond_list_a2,
                        const fp_type* const __restrict__* const __restrict__ grid_maps,
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
    constexpr_for<0, get_num_atom_clusters(), 1>([&](const auto cluster_index) {
      const auto num_atoms_cluster_prev = atoms_clusters[cluster_index - 1];
      const auto num_atoms_cluster      = atoms_clusters[cluster_index];
      if (num_atoms < num_atoms_cluster && num_atoms >= num_atoms_cluster_prev)
        // Simulate the population evolution for the given amount of time
        evaluate_fitness_impl<num_atoms_cluster>(ligand_x,
                                                 ligand_y,
                                                 ligand_z,
                                                 ligand_vol,
                                                 ligand_solpar,
                                                 ligand_charge,
                                                 ligand_num_hbond,
                                                 ligand_Rij_hb,
                                                 ligand_Rii,
                                                 ligand_epsij_hb,
                                                 ligand_epsii,
                                                 map_ligand_types,
                                                 num_atoms,
                                                 num_rotamers,
                                                 frag_masks,
                                                 frag_start_indexes,
                                                 frag_stop_indexes,
                                                 num_nonbond,
                                                 non_bond_list_a1,
                                                 non_bond_list_a2,
                                                 grid_maps,
                                                 electro_map,
                                                 desolv_map,
                                                 num_generations,
                                                 population_size,
                                                 tournament_length,
                                                 mutation_prob,
                                                 minimum,
                                                 maximum,
                                                 center,
                                                 map_index_x,
                                                 map_index_xy,
                                                 population_buffer1,
                                                 population_buffer2,
                                                 seed);
    });
  }
} // namespace mudock
