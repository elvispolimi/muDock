#pragma once

#include <memory>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/grid.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <random>
#include <span>
#include <string>
#include <vector>

namespace mudock {

  /**
   * The virtual screening algorithm is basically a genetic algorithm that use the ligand
   * energy as fitness function, and geometric transformations of the molecule as genes.
   */
  class virtual_screen_cpp {
    // these are information about the target protein
    std::shared_ptr<const grid_atom_mapper> grid_atom_maps;
    std::shared_ptr<const grid_map> electro_map;
    std::shared_ptr<const grid_map> desolv_map;

    // define the GA population
    std::vector<individual> population;
    std::vector<individual> next_population;

    // this algorithm requires a random source to work
    std::mt19937 generator;
    std::uniform_real_distribution<fp_type> dist{fp_type{0.0}, fp_type{1.0}};

    // TODO
    // random_generator<int> rnd_gen;

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
    virtual_screen_cpp(std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                       std::shared_ptr<const grid_map>& electro_map,
                       std::shared_ptr<const grid_map>& desolv_map,
                       const knobs& knobs);

    void operator()(static_molecule& ligand);
  };

} // namespace mudock
