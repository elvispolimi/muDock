#include "mudock/grid/point3D.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <mudock/chem.hpp>
#include <mudock/grid.hpp>
#include <mudock/type_alias.hpp>
#include <vector>

namespace mudock {
  fp_type calc_ddd_Mehler_Solmajer(fp_type distance) {
    /*____________________________________________________________________________
     * Distance-dependent dielectric ewds: Mehler and Solmajer, Prot Eng 4, 903-910.
     *____________________________________________________________________________*/
    const fp_type lambda{0.003627};
    const fp_type epsilon0{78.4};
    const fp_type A{-8.5525};
    const fp_type B = epsilon0 - A;
    const fp_type rk{7.7839};
    const fp_type lambda_B = -lambda * B;

    fp_type epsilon = A + B / (fp_type{1} + rk * std::exp(lambda_B * distance));

    if (epsilon < std::numeric_limits<fp_type>::epsilon()) {
      epsilon = 1.0L;
    }
    return epsilon;
  }

  grid_map generate_electrostatic_grid_map(dynamic_molecule& receptor) {
    //  Get maximum and minimum of the bounding box around the receptor
    const fp_type receptor_max_x = std::ranges::max(receptor.get_x());
    const fp_type receptor_max_y = std::ranges::max(receptor.get_y());
    const fp_type receptor_max_z = std::ranges::max(receptor.get_z());
    const fp_type receptor_min_x = std::ranges::min(receptor.get_x());
    const fp_type receptor_min_y = std::ranges::min(receptor.get_y());
    const fp_type receptor_min_z = std::ranges::min(receptor.get_z());

    const index3D npts{
        static_cast<size_t>(ceilf((receptor_max_x - receptor_min_x + cutoff_distance * 2) / grid_spacing)),
        static_cast<size_t>(ceilf((receptor_max_y - receptor_min_y + cutoff_distance * 2) / grid_spacing)),
        static_cast<size_t>(ceilf((receptor_max_z - receptor_min_z + cutoff_distance * 2) / grid_spacing))};
    const point3D grid_minimum{floorf(receptor_min_x - cutoff_distance),
                               floorf(receptor_min_y - cutoff_distance),
                               floorf(receptor_min_z - cutoff_distance)};

    grid_map electrostatic_map{npts};

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
