#pragma once

#include "mudock/grid/grid_map.hpp"

#include <mudock/grid.hpp>
#include <mudock/molecule.hpp>
#include <vector>

namespace mudock {

  void calc_energy(const static_molecule& receptor,
                   const dynamic_molecule& ligand,
                   const grid_atom_mapper& grid_atom_maps,
                   const grid_map& electro_map,
                   const grid_map& desolv_map);

}
