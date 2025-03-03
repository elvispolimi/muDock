#include "mudock/grid/mdindex.hpp"

#include <algorithm>
#include <array>
#include <mudock/chem.hpp>
#include <mudock/grid.hpp>
#include <vector>

namespace mudock {
  grid_map generate_desolvation_grid_map(dynamic_molecule& receptor) {
    grid_map desolvation_map{receptor};
    // Get grid infos from the first element
    const auto& grid_minimum = desolvation_map.minimum;
    const auto& npts         = desolvation_map.index;

    /* exponential function for receptor and ligand desolvation */
    /* note: the solvation term ranges beyond the non-bond cutoff 
    * and will not be smoothed 
    */
    std::array<fp_type, NDIEL> sol_fn;
    for (int indx_r = 1; indx_r < NDIEL; indx_r++) {
      const fp_type r = indx_r / A_DIV;
      sol_fn[indx_r] =
          autodock_parameters::coeff_desolv * std::exp(-(r * r) / (fp_type{2} * (sigma * sigma)));
    }

    for (int index_z = 0; index_z < npts.size_z(); ++index_z) {
      /*
      *  c[0:2] contains the current grid point.
      */
      const fp_type coord_z = grid_minimum.z + index_z * grid_spacing;
      for (int index_y = 0; index_y < npts.size_y(); ++index_y) {
        const fp_type coord_y = grid_minimum.y + index_y * grid_spacing;
        for (int index_x = 0; index_x < npts.size_x(); ++index_x) {
          const fp_type coord_x = grid_minimum.x + index_x * grid_spacing;
          fp_type energy{0};
          for (int index = 0; index < receptor.num_atoms(); ++index) {
            const fp_type d = distance(point3D{receptor.x(index), receptor.y(index), receptor.z(index)},
                                       point3D{coord_x, coord_y, coord_z});
            if (d > NBC)
              continue; /* onto the next atom... */
            const int indx_r = std::min<int>(std::floor(d * A_DIV), NDIEL - 1);
            energy += solpar_q * receptor.get_vol()[index] * sol_fn[indx_r];
          }
          desolvation_map.at(index_x, index_y, index_z) = energy;
        }
      }
    }

    return desolvation_map;
  }

} // namespace mudock
