#pragma once

#include <mudock/grid.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/fragments.hpp>
#include <vector>

namespace mudock {

  fp_type calc_energy(const dynamic_molecule& receptor,
                      const static_molecule& ligand,
                      const fragments<static_containers>& ligand_fragments,
                      const grid_atom_mapper& grid_atom_maps,
                      const grid_map& electro_map,
                      const grid_map& desolv_map);
}
