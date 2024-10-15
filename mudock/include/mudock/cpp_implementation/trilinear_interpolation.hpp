#pragma once

#include <mudock/grid.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  inline fp_type trilinear_interpolation(const grid_map& map, const point3D& coord) {
    const point3D coordinates_in_grid = map.get_index_from_coordinates(coord);
    const point3D coordinates_floored = point3D{std::floor(coordinates_in_grid.x),
                                                std::floor(coordinates_in_grid.y),
                                                std::floor(coordinates_in_grid.z)};

    const fp_type delta_low_x = coordinates_in_grid.x - coordinates_floored.x;
    const fp_type delta_low_y = coordinates_in_grid.y - coordinates_floored.y;
    const fp_type delta_low_z = coordinates_in_grid.z - coordinates_floored.z;

    std::array<fp_type, 2> px{fp_type{1} - delta_low_x, delta_low_x};
    std::array<fp_type, 2> py{fp_type{1} - delta_low_y, delta_low_y};
    std::array<fp_type, 2> pz{fp_type{1} - delta_low_z, delta_low_z};
    fp_type value{0};
    for (int i = 0; i <= 1; ++i)
      for (int j = 0; j <= 1; ++j)
        for (int t = 0; t <= 1; ++t) {
          value += px[t] * py[j] * pz[i] *
                   map.at(coordinates_floored.x + t, coordinates_floored.y + j, coordinates_floored.z + i);
        }
    return value;
  }
} // namespace mudock
