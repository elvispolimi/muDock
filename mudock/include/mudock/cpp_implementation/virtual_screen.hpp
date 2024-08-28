#pragma once

#include <memory>
#include <mudock/grid.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <random>
#include <span>
#include <string>

namespace mudock {

  class virtual_screen_cpp {
    // this will be the scratchpad memory for the cpp implementation
    // NOTE: we will implement a random scoring function to test the infrastructure
    std::mt19937 generator;
    std::uniform_real_distribution<fp_type> dist;

    std::shared_ptr<const dynamic_molecule>& protein;
    std::shared_ptr<const grid_atom_mapper> grid_atom_maps;
    std::shared_ptr<const grid_map> electro_map;
    std::shared_ptr<const grid_map> desolv_map;


    // the configuration of the GA algorithm
    knobs configuration;

  public:
    virtual_screen_cpp(std::shared_ptr<const dynamic_molecule>& protein,
                       std::shared_ptr<const grid_atom_mapper> grid_atom_maps,
                       std::shared_ptr<const grid_map> electro_map,
                       std::shared_ptr<const grid_map> desolv_map,
                       const knobs knobs);

    void operator()(static_molecule& ligand);
  };

} // namespace mudock
