#include "mudock/grid/mdindex.hpp"

#include <algorithm>
#include <array>
#include <mudock/chem.hpp>
#include <mudock/grid.hpp>
#include <vector>

namespace mudock {
  grid_map generate_desolvation_grid_map(dynamic_molecule& receptor) {
    //  Get maximum and minimum of the bounding box around the receptor
    const fp_type max_x = std::ranges::max(receptor.get_x());
    const fp_type max_y = std::ranges::max(receptor.get_y());
    const fp_type max_z = std::ranges::max(receptor.get_z());
    const fp_type min_x = std::ranges::min(receptor.get_x());
    const fp_type min_y = std::ranges::min(receptor.get_y());
    const fp_type min_z = std::ranges::min(receptor.get_z());

    const point3D p_center{(max_x - min_x) / fp_type{2},
                           (max_y - min_y) / fp_type{2},
                           (max_z - min_z) / fp_type{2}};
    const index3D npts{static_cast<size_t>((max_x - min_x) / grid_spacing + 1),
                       static_cast<size_t>((max_y - min_y) / grid_spacing + 1),
                       static_cast<size_t>((max_z - min_z) / grid_spacing + 1)};

    grid_map desolvation_map{npts};

    // TODO

    return desolvation_map;
  }

} // namespace mudock
