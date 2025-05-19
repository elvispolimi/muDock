#pragma once

#include <mudock/chem/autodock_ligand_types.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/grid.hpp>
#include <mudock/grid/mdspan.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  static constexpr auto get_num_map = num_autodock_ligand_types() + 2;

  /**
   * This data structure represents how autodock encode the target protein. The idea is to represent different
   * the 3D space using a set discretized grids. Each grid model one property. We always have two grid to
   * represent electrostatic and desolvation components. Then, we can have a set of grid that represent a
   * contribution for different ligand's atom types.
   */
  struct autodock_protein {
    inline auto& get_eletrostatic() { return map[0]; };
    inline auto& get_desolvation() { return map[1]; };
    inline auto& get_atom_map(const autodock_ligand_ff type) { return map[static_cast<int>(type)]; };

    autodock_protein(const point3D min, const point3D max, const fp_type resolution);

  private:
    md_index<3> index;
    md_container<std::vector<fp_type>, 1> data;
    std::array<space_grid_view, get_num_map> map;
  };

  // this function will generate the autodock grids from the parsed protein
  autodock_protein make_autodock_protein(const dynamic_molecule& protein, const molecule_graph_type& graph);

} // namespace mudock
