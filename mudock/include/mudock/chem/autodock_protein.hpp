#pragma once

#include <initializer_list>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/grid.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <unordered_map>

namespace mudock {

  /**
   * This data structure represents how autodock encode the target protein. The idea is to represent different
   * the 3D space using a set discretized grids. Each grid model one property. We always have two grid to
   * represent electrostatic and desolvation components. Then, we can have a set of grid that represent a
   * contribution for different ligand's atom types.
   */
  struct autodock_protein {
    space_grid electrostatic;
    space_grid desolvation;
    std::unordered_map<autodock_ff, space_grid> atom_map;
  };

  // this function will generate the autodock grids from the parsed protein
  autodock_protein make_autodock_protein(const dynamic_molecule& protein, const molecule_graph_type& graph);

} // namespace mudock
