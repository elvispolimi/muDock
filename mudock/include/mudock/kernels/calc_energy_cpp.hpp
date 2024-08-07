#pragma once

#include <mudock/grid.hpp>
#include <mudock/molecule.hpp>
#include <vector>

namespace mudock {

  void calc_energy(const static_molecule& receptor, const dynamic_molecule& ligand, const std::vector<grid_map> grid_maps);

}
