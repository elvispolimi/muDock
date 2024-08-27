#include <algorithm>
#include <array>
#include <limits>
#include <mudock/chem.hpp>
#include <mudock/grid.hpp>
#include <mudock/type_alias.hpp>
#include <vector>

namespace mudock {
  grid_map generate_electrostatic_grid_map(dynamic_molecule& receptor) {
    //  Get maximum and minimum of the bounding box around the receptor
    const fp_type receptor_max_x = std::ranges::max(receptor.get_x());
    const fp_type receptor_max_y = std::ranges::max(receptor.get_y());
    const fp_type receptor_max_z = std::ranges::max(receptor.get_z());
    const fp_type receptor_min_x = std::ranges::min(receptor.get_x());
    const fp_type receptor_min_y = std::ranges::min(receptor.get_y());
    const fp_type receptor_min_z = std::ranges::min(receptor.get_z());

    const point3D grid_minimum{std::floor((receptor_min_x - cutoff_distance) * fp_type{2}) / fp_type{2},
                               std::floor((receptor_min_y - cutoff_distance) * fp_type{2}) / fp_type{2},
                               std::floor((receptor_min_z - cutoff_distance) * fp_type{2}) / fp_type{2}};
    const point3D grid_maximum{std::ceil((receptor_max_x + cutoff_distance) * fp_type{2}) / fp_type{2},
                               std::ceil((receptor_max_y + cutoff_distance) * fp_type{2}) / fp_type{2},
                               std::ceil((receptor_max_z + cutoff_distance) * fp_type{2}) / fp_type{2}};
    const index3D npts{static_cast<size_t>((grid_maximum.x - grid_minimum.x) / grid_spacing),
                       static_cast<size_t>((grid_maximum.y - grid_minimum.y) / grid_spacing),
                       static_cast<size_t>((grid_maximum.z - grid_minimum.z) / grid_spacing)};

    grid_map electrostatic_map{npts, grid_minimum, grid_maximum};

    std::array<fp_type, NDIEL> epsilon_fn;
    std::array<fp_type, NDIEL> r_epsilon_fn;
    epsilon_fn[0] = 1.0;
    for (size_t indx_r = 1; indx_r < NDIEL; ++indx_r) {
      epsilon_fn[indx_r] = calc_ddd_Mehler_Solmajer(indx_r / A_DIV);
    }
    /* convert epsilon to factor / epsilon */
    for (size_t i = 0; i < NDIEL; ++i) { r_epsilon_fn[i] = factor / epsilon_fn[i]; }

    for (size_t index_z = 0; index_z < npts.size_z(); ++index_z) {
      /*
      *  c[0:2] contains the current grid point.
      */
      const fp_type coord_z = grid_minimum.z + index_z * grid_spacing;
      for (size_t index_y = 0; index_y < npts.size_y(); ++index_y) {
        const fp_type coord_y = grid_minimum.y + index_y * grid_spacing;
        for (size_t index_x = 0; index_x < npts.size_x(); ++index_x) {
          const fp_type coord_x = grid_minimum.x + index_x * grid_spacing;
          fp_type energy{0};
          for (size_t index = 0; index < receptor.num_atoms(); ++index) {
            const fp_type d = distance(point3D{receptor.x(index), receptor.y(index), receptor.z(index)},
                                       point3D{coord_x, coord_y, coord_z});
            const fp_type inv_rmax = fp_type{1} / std::max(d, fp_type{0.5});
            const size_t indx_r    = std::min<size_t>(std::floor(d * A_DIV), NDIEL - 1);
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
