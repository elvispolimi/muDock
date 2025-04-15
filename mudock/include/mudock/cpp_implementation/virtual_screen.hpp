#pragma once

#include <memory>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cpp_implementation/evaluate_fitness_cpp.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/grid.hpp>
#include <mudock/grid/grid_map.hpp>
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
    std::shared_ptr<const grid_atom_mapper> grid_atom_maps;
    std::shared_ptr<const grid_map> electro_map;
    std::shared_ptr<const grid_map> desolv_map;

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
    virtual_screen_cpp(std::shared_ptr<const grid_atom_mapper>& _grid_atom_maps,
                       std::shared_ptr<const grid_map>& _electro_map,
                       std::shared_ptr<const grid_map>& _desolv_map,
                       const knobs& knobs)
        : grid_atom_maps(_grid_atom_maps),
          electro_map(_electro_map),
          desolv_map(_desolv_map),
          population(knobs.population_number),
          next_population(knobs.population_number),
          configuration(knobs) {}

    void operator()(static_molecule& ligand) {
      const auto seed =
          configuration.seed.has_value()
              ? configuration.seed.value()
              : static_cast<size_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
      // Place the molecule to the center of the target protein
      const int num_atoms = ligand.num_atoms();
      const auto x = ligand.get_x(), y = ligand.get_y(), z = ligand.get_z();
      const auto ligand_center_of_mass = compute_center_of_mass(x, y, z);
      translate_molecule<cpu_vectorization::AUTO>(x.data(),
                                                  y.data(),
                                                  z.data(),
                                                  num_atoms,
                                                  electro_map->center.x - ligand_center_of_mass.x,
                                                  electro_map->center.y - ligand_center_of_mass.y,
                                                  electro_map->center.z - ligand_center_of_mass.z);

      const auto num_rotamers = ligand.num_rotamers();

      // Get weed bonds and non bonds lists
      std::vector<int> non_bond_list_a1, non_bond_list_a2;
      std::vector<mudock::fp_type> cA_v, cB_v;
      std::vector<int> xB_v;
      mudock::non_bond_list(ligand, non_bond_list_a1, non_bond_list_a2);
      mudock::precompute_lennard_jones(non_bond_list_a1.size(),
                                       cA_v,
                                       cB_v,
                                       xB_v,
                                       ligand,
                                       non_bond_list_a1,
                                       non_bond_list_a2);

      const int atom_map_size         = grid_atom_maps.get()->get_single_map_size();
      const fp_type* atom_map_pointer = grid_atom_maps.get()->get_fused_maps().data();

      std::vector<int> frag_masks;
      std::vector<int> frag_start_indexes;
      std::vector<int> frag_stop_indexes;
      get_linearized_fragments_mask(num_atoms,
                                    num_rotamers,
                                    frag_masks,
                                    frag_start_indexes,
                                    frag_stop_indexes,
                                    ligand);

      std::vector<int> map_ligand_offsets;
      map_ligand_offsets.resize(num_atoms);
      for (int i = 0; i < num_atoms; i++)
        map_ligand_offsets[i] =
            static_cast<int>(map_from_autodock_type(ligand.autodock_type(i))) * atom_map_size;
      // Simulate the population evolution for the given amount of time
      evaluate_fitness<vect>(x.data(),
                             y.data(),
                             z.data(),
                             ligand.get_vol().data(),
                             ligand.get_solpar().data(),
                             ligand.get_charge().data(),
                             map_ligand_offsets.data(),
                             num_atoms,
                             num_rotamers,
                             frag_masks.data(),
                             frag_start_indexes.data(),
                             frag_stop_indexes.data(),
                             non_bond_list_a1.size(),
                             non_bond_list_a1.data(),
                             non_bond_list_a2.data(),
                             cA_v.data(),
                             cB_v.data(),
                             xB_v.data(),
                             atom_map_pointer,
                             electro_map.get()->data(),
                             desolv_map.get()->data(),
                             configuration.num_generations,
                             configuration.population_number,
                             configuration.tournament_length,
                             configuration.mutation_prob,
                             electro_map.get()->minimum_coord.get_array().data(),
                             electro_map.get()->maximum_coord.get_array().data(),
                             electro_map.get()->center.get_array().data(),
                             electro_map.get()->index.size_x(),
                             electro_map.get()->index.size_xy(),
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
                                     frag_masks.data(),
                                     frag_start_indexes.data(),
                                     frag_stop_indexes.data());
      ligand.properties.assign(property_type::SCORE, std::to_string(best_individual_it->score));
    }
  };

} // namespace mudock
