#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/ligand_maps.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/cpp_implementation/calc_energy_cpp.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/trilinear_interpolation.hpp>
#include <mudock/utils.hpp>
#include <random>
#include <stdexcept>
#include <mudock/fjapp_utils.hpp>

#define FLATTENED_2D(x, y, index_x) ((y) * index_x + (x))

namespace mudock {
  static constexpr auto coordinate_step = fp_type{0.2};
  static constexpr auto angle_step      = fp_type{4};

  template<typename T>
  [[nodiscard]] const T random_gen_cpp(std::mt19937& generator,
                                       std::uniform_real_distribution<fp_type>& dist,
                                       const T& min,
                                       const T& max) {
    fp_type value;
    if constexpr (is_debug())
      // TODO value here for debug
      value = fp_type{0.4};
    else {
      value = dist(generator);
    }
    return static_cast<T>(value * (max - min) + min);
  }

  int get_selection_distribution(std::mt19937& generator,
                                 std::uniform_real_distribution<fp_type>& dist,
                                 const int& population_number) {
    return random_gen_cpp<int>(generator, dist, 0, population_number - 1);
  };
  fp_type get_init_change_distribution(std::mt19937& generator,
                                       std::uniform_real_distribution<fp_type>& dist) {
    return random_gen_cpp<fp_type>(generator, dist, -45, 45);
  }
  fp_type get_mutation_change_distribution(std::mt19937& generator,
                                           std::uniform_real_distribution<fp_type>& dist) {
    return random_gen_cpp<fp_type>(generator, dist, -10, 10);
  };
  fp_type get_mutation_coin_distribution(std::mt19937& generator,
                                         std::uniform_real_distribution<fp_type>& dist) {
    return random_gen_cpp<fp_type>(generator, dist, 0, 1);
  };
  int get_crossover_distribution(std::mt19937& generator,
                                 std::uniform_real_distribution<fp_type>& dist,
                                 const int& num_rotamers) {
    return random_gen_cpp<int>(generator, dist, 0, 6 + num_rotamers);
  };
  const chromosome& tournament_selection(std::mt19937& generator,
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

  fp_type calc_energy(const fp_type* __restrict__ ligand_x,
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
                      const autodock_ff* __restrict__ ligand_autodock_type,
                      const int num_atoms,
                      const int n_torsions,
                      const uint_fast8_t* __restrict__ nbmatrix,
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
    for (int index = 0; index < num_atoms; ++index) {
      fp_type coord[3]{ligand_x[index], ligand_y[index], ligand_z[index]};
      const auto& atom_charge = ligand_charge[index];

      if (coord[0] < minimum[0] || coord[0] > maximum[0] || coord[1] < minimum[1] || coord[1] > maximum[1] ||
          coord[2] < minimum[2] || coord[2] > maximum[2]) {
        // printf("Atom %d is outside\n", index);
        const fp_type dist = std::pow(std::fabs(coord[0] - center[0]), fp_type{2}) +
                             std::pow(std::fabs(coord[1] - center[1]), fp_type{2}) +
                             std::pow(std::fabs(coord[2] - center[2]), fp_type{2});
        const fp_type epenalty = dist * ENERGYPENALTY;
        elect_total_trilinear += epenalty;
        emap_total_trilinear += epenalty;
      } else {
        const fp_type* atom_map =
            grid_maps[static_cast<int>(map_from_autodock_type(ligand_autodock_type[index]))];

        coord[0] = (coord[0] - minimum[0]) * inv_spacing;
        coord[1] = (coord[1] - minimum[1]) * inv_spacing;
        coord[2] = (coord[2] - minimum[2]) * inv_spacing;

        // Trilinear Interpolationp
        elect_total_trilinear +=
            trilinear_interpolation(electro_map, coord, map_index_x, map_index_xy) * atom_charge;
        emap_total_trilinear += trilinear_interpolation(atom_map, coord, map_index_x, map_index_xy);
        dmap_total_trilinear +=
            trilinear_interpolation(desolv_map, coord, map_index_x, map_index_xy) * std::fabs(atom_charge);
      }
    }

    fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
    if (n_torsions > 0) {
      // TODO @Davide suppose that the receptor does not have Flexible residues eintcal.cc:147
      // TODO

      for (int i = 0; i < num_atoms; ++i)
        for (int j = i + 1; j < num_atoms; ++j) {
          // for (const auto& non_bond: non_bond_list) {
          if ((nbmatrix[FLATTENED_2D(i, j, num_atoms)] == 1 &&
               nbmatrix[FLATTENED_2D(j, i, num_atoms)] == 1)) {
            const int& a1 = i;
            const int& a2 = j;

            const fp_type distance_two = std::pow(std::fabs(ligand_x[a1] - ligand_x[a2]), fp_type{2}) +
                                         std::pow(std::fabs(ligand_y[a1] - ligand_y[a2]), fp_type{2}) +
                                         std::pow(std::fabs(ligand_z[a1] - ligand_z[a2]), fp_type{2});
            const fp_type distance_two_clamp = std::clamp(distance_two, RMIN_ELEC * RMIN_ELEC, distance_two);
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
                                     std::exp(fp_type{-0.5} / (sigma * sigma) * distance_two_clamp) *
                                     nb_desolv;
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

                /* smooth with min function; 
              r_smooth is Angstrom range of "smoothing" */
                // TODO rework considering this precomputing with other values
                // if (r_smooth > 0) {
                //   const fp_type rlow  = distance - r_smooth / 2;
                //   const fp_type rhigh = distance + r_smooth / 2;
                //   fp_type energy_smooth{100000};
                //   // ((((int)((r)*A_DIV)) > NEINT_1) ? NEINT_1 : ((int)((r)*A_DIV))
                //   for (int j = std::max(fp_type{0}, std::min(rlow * A_DIV, fp_type{NEINT - 1}));
                //        j <= std::min(fp_type{NEINT - 1}, std::min(rhigh * A_DIV, fp_type{NEINT - 1}));
                //        ++j)
                //     energy_smooth = std::min(energy_smooth, e_vdW_Hb);

                //   e_vdW_Hb = energy_smooth;
                // } /* endif smoothing */

                // TODO check energy smoothing intnbtable.cc:215
                // with NOSQRT to False it seems to be the same as calmping to EINTCLAMP
              } else {
                throw std::runtime_error("ERROR: Exponents must be different, to avoid division by zero!");
              }
            }
            emap_total_eintcal += e_vdW_Hb;
          }
          // else if ((nbmatrix.at(i, j) != 0 && nbmatrix.at(j, i) == 0) ||
          //                (nbmatrix.at(i, j) == 0 && nbmatrix.at(j, i) != 0)) {
          //       std::ostringstream oss;
          //       // Build the formatted string
          //       oss << "BUG: ASSYMMETRY detected in Non-Bond Matrix at " << i << "," << j;
          //       error(oss.str());
          //     }
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
                        const autodock_ff* __restrict__ ligand_autodock_type,
                        const int num_atoms,
                        const int num_rotamers,
                        const int* __restrict__ frag_masks,
                        const int* __restrict__ frag_start_indexes,
                        const int* __restrict__ frag_stop_indexes,
                        const uint_fast8_t* __restrict__ nbmatrix,
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

    FJAPP_MARKER_START("GA");

    auto* population      = population_buffer1;
    auto* next_population = population_buffer2;
    auto altered_x        = std::make_unique<std::array<fp_type, max_static_atoms()>>();
    auto altered_y        = std::make_unique<std::array<fp_type, max_static_atoms()>>();
    auto altered_z        = std::make_unique<std::array<fp_type, max_static_atoms()>>();

    // Randomly initialize the population
    for (int element_index = 0; element_index < population_size; ++element_index) {
      auto& element = population[element_index];
      for (int i{0}; i < 3; ++i) { // initialize the rigid translation
        element.genes[i] = get_init_change_distribution(generator, dist) * coordinate_step;
      }
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
                                        ligand_autodock_type,
                                        num_atoms,
                                        num_rotamers,
                                        nbmatrix,
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
  }
} // namespace mudock
