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
    std::shared_ptr<const dynamic_molecule>& protein;
    std::shared_ptr<const grid_atom_mapper> grid_atom_maps;
    std::shared_ptr<const grid_map> electro_map;
    std::shared_ptr<const grid_map> desolv_map;

    // define the GA population
    struct individual {
      chromosome genes;
      fp_type score;
    };
    std::vector<individual> population;
    std::vector<individual> next_population;

    // this algorithm requires a random source to work
    std::mt19937 generator;

    // the configuration of the GA algorithm
    knobs configuration;

    // utility function to select a parent for the crossover
    const chromosome& tournament_selection();

  public:
    virtual_screen_cpp(std::shared_ptr<const dynamic_molecule>& protein,
                       std::shared_ptr<const grid_atom_mapper> grid_atom_maps,
                       std::shared_ptr<const grid_map> electro_map,
                       std::shared_ptr<const grid_map> desolv_map,
                       const knobs knobs);

    void operator()(static_molecule& ligand);
  };

} // namespace mudock
