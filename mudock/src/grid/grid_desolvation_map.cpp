#include "mudock/grid/mdindex.hpp"

#include <algorithm>
#include <array>
#include <mudock/chem.hpp>
#include <mudock/grid.hpp>
#include <vector>

namespace mudock {
  grid_map generate_desolvation_grid_map(dynamic_molecule& receptor) {
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

    grid_map desolvation_map{npts};

    /* exponential function for receptor and ligand desolvation */
    /* note: the solvation term ranges beyond the non-bond cutoff 
    * and will not be smoothed 
    */
    std::array<fp_type, NDIEL> sol_fn;
    for (size_t indx_r = 1; indx_r < NDIEL; indx_r++) {
      const fp_type r = indx_r / A_DIV;
      sol_fn[indx_r] =
          autodock_parameters::coeff_desolv * std::exp(-std::sqrt(r) / (fp_type{2} * std::sqrt(sigma)));
    }

    for (size_t index_z = 0; index_z < npts.size_z(); ++index_z) {
      /*
      *  c[0:2] contains the current grid point.
      */
      const fp_type coord_z = grid_minimum.z + index_z * grid_spacing;
      for (size_t index_y = 0; index_y < npts.size_z(); ++index_y) {
        const fp_type coord_y = grid_minimum.y + index_y * grid_spacing;
        for (size_t index_x = 0; index_x < npts.size_z(); ++index_x) {
          const fp_type coord_x = grid_minimum.x + index_x * grid_spacing;
          fp_type energy{0};
          for (size_t index = 0; index < receptor.num_atoms(); ++index) {
            point3D dist = difference(point3D{receptor.x(index), receptor.y(index), receptor.z(index)},
                                      point3D{coord_x, coord_y, coord_z});
            fp_type d    = sqrtf(sum_components(square(dist)));
            //  TODO check @Davide why no error here?
            if (d == fp_type{0}) {
              d = std::numeric_limits<fp_type>::epsilon();
            }
            const size_t indx_r = std::min<size_t>(std::floor(d * A_DIV), NDIEL - 1);

            energy += solpar_q * receptor.get_vol()[index] * sol_fn[indx_r];
          }
          desolvation_map.at(coord_x, coord_y, coord_z) = energy;
        }
      }
    }

    return desolvation_map;
  }

} // namespace mudock
