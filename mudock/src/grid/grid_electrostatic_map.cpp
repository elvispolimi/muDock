#include <algorithm>
#include <array>
#include <limits>
#include <mudock/chem.hpp>
#include <mudock/grid.hpp>
#include <mudock/type_alias.hpp>
#include <vector>

namespace mudock {
  grid_map generate_electrostatic_grid_map(dynamic_molecule& receptor) {
    grid_map electrostatic_map{receptor};
    // Get grid infos from the first element
    const auto& grid_minimum = electrostatic_map.minimum;
    const auto& npts         = electrostatic_map.index;

    std::array<fp_type, NDIEL> epsilon_fn;
    std::array<fp_type, NDIEL> r_epsilon_fn;
    epsilon_fn[0] = 1.0;
    for (int indx_r = 1; indx_r < NDIEL; ++indx_r) {
      epsilon_fn[indx_r] = calc_ddd_Mehler_Solmajer(indx_r / A_DIV);
    }
    /* convert epsilon to factor / epsilon */
    for (int i = 0; i < NDIEL; ++i) { r_epsilon_fn[i] = factor / epsilon_fn[i]; }

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
            const fp_type inv_rmax = fp_type{1} / std::max(d, fp_type{0.5});
            const int indx_r       = std::min<int>(std::floor(d * A_DIV), NDIEL - 1);
            energy +=
                receptor.charge(index) * inv_rmax * r_epsilon_fn[indx_r] * autodock_parameters::coeff_estat;
          }
          electrostatic_map.at(index_x, index_y, index_z) = energy;
        }
      }
    }

    return electrostatic_map;
  }

} // namespace mudock
