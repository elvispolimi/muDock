#pragma once

#include "mudock/chem/autodock_ligand.hpp"

#include <mudock/chem/autodock_protein.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cpp_implementation/evaluate_fitness_cpp.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/grid.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/fragments.hpp>
#include <mudock/type_alias.hpp>
#include <vector>

namespace mudock {
  /**
   * The virtual screening algorithm is basically a genetic algorithm that use the ligand
   * energy as fitness function, and geometric transformations of the molecule as genes.
   */
  template<cpu_vectorization vect>
  class virtual_screen_cpp {
    // these are information about the target protein
    const autodock_protein& adt_protein;

    // define the GA population
    std::vector<individual> population;
    std::vector<individual> next_population;

    // the configuration of the GA algorithm
    knobs configuration;

    // utility function to select a parent for the crossover
    const chromosome& tournament_selection();

    template<typename T>
    [[nodiscard]] const T random_gen_cpp(const T& min, const T& max);

    [[nodiscard]] int get_selection_distribution();
    [[nodiscard]] fp_type get_init_change_distribution();
    [[nodiscard]] fp_type get_mutation_change_distribution();
    [[nodiscard]] fp_type get_mutation_coin_distribution();
    [[nodiscard]] int get_crossover_distribution(const int& num_rotamers);

  public:
    virtual_screen_cpp(const autodock_protein& _adt_protein, const knobs& knobs)
        : adt_protein(_adt_protein),
          population(knobs.population_number),
          next_population(knobs.population_number),
          configuration(knobs) {}

    void operator()(static_molecule& ligand) {
      const auto seed =
          configuration.seed.has_value()
              ? configuration.seed.value()
              : static_cast<size_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
      // Place the molecule to the center of the target protein
      const int num_atoms     = ligand.num_atoms();
      const auto num_rotamers = ligand.num_rotamers();

      auto x = ligand.get_x(), y = ligand.get_y(), z = ligand.get_z();
      const auto ligand_center_of_mass = compute_center_of_mass(x, y, z);
      const auto offset                = adt_protein.get_center() - ligand_center_of_mass;
      translate_molecule<cpu_vectorization::AUTO>(x.data(),
                                                  y.data(),
                                                  z.data(),
                                                  num_atoms,
                                                  offset.x(),
                                                  offset.y(),
                                                  offset.z());
      // for (int a = 0; a < num_atoms; ++a) printf("%f %f %f\n", x[a], y[a], z[a]);

      auto adt_ligand = autodock_ligand{ligand};
      adt_ligand.update_offsets(adt_protein);

      // Simulate the population evolution for the given amount of time
      evaluate_fitness<vect>(adt_ligand.get_ligand_x_p(),
                             adt_ligand.get_ligand_y_p(),
                             adt_ligand.get_ligand_z_p(),
                             adt_ligand.get_ligand_vol(),
                             adt_ligand.get_ligand_solpar(),
                             adt_ligand.get_ligand_charge(),
                             adt_ligand.get_atom_map_offsets(),
                             num_atoms,
                             num_rotamers,
                             adt_ligand.get_fragments_masks(),
                             adt_ligand.get_fragmets_starts(),
                             adt_ligand.get_fragments_stops(),
                             adt_ligand.get_non_bond_size(),
                             adt_ligand.get_non_bond_A(),
                             adt_ligand.get_non_bond_B(),
                             adt_ligand.get_non_bond_cA(),
                             adt_ligand.get_non_bond_cB(),
                             adt_ligand.get_non_bond_xB(),
                             adt_protein.get_maps_pointer(),
                             configuration.num_generations,
                             configuration.population_number,
                             configuration.tournament_length,
                             configuration.mutation_prob,
                             adt_protein.get_min_p(),
                             adt_protein.get_max_p(),
                             adt_protein.get_center_p(),
                             adt_protein.get_size_x(),
                             adt_protein.get_size_xy(),
                             adt_protein.get_size_xyz(),
                             population.data(),
                             next_population.data(),
                             seed);

      // update the ligand position with the best one that we found
      const std::vector<individual>& last_population =
          (configuration.num_generations % 2 == 0) ? next_population : population;

      const auto best_individual_it =
          std::min_element(std::begin(last_population),
                           std::end(last_population),
                           [](const auto a, const auto b) { return a.score < b.score; });
      apply<cpu_vectorization::AUTO>(x.data(),
                                     y.data(),
                                     z.data(),
                                     best_individual_it->genes,
                                     num_atoms,
                                     num_rotamers,
                                     adt_ligand.get_fragments_masks(),
                                     adt_ligand.get_fragmets_starts(),
                                     adt_ligand.get_fragments_stops());
      ligand.properties.assign(property_type::SCORE, std::to_string(best_individual_it->score));
    }
  };

} // namespace mudock
