#pragma once

#include <mudock/type_alias.hpp>

#define FLATTENED_3D(x, y, z, index_x, index_xy) (index_xy * (z) + (y) * index_x + (x))

namespace mudock {
  inline fp_type trilinear_interpolation(const fp_type* __restrict__ map,
                                         const fp_type* __restrict__ coord,
                                         const int& map_index_x,
                                         const int& map_index_xy) {
    const int u0      = coord[0];
    const fp_type p0u = coord[0] - static_cast<fp_type>(u0);
    const fp_type p1u = fp_type{1} - p0u;

    const int v0      = coord[1];
    const fp_type p0v = coord[1] - static_cast<fp_type>(v0);
    const fp_type p1v = fp_type{1} - p0v;

    const int w0      = coord[2];
    const fp_type p0w = coord[2] - static_cast<fp_type>(w0);
    const fp_type p1w = fp_type{1} - p0w;

    const fp_type pu[2] = {p1u, p0u};
    const fp_type pv[2] = {p1v, p0v};
    const fp_type pw[2] = {p1w, p0w};
    fp_type value{0};
#pragma unroll
    for (int i = 0; i <= 1; i++)
#pragma unroll
      for (int t = 0; t <= 1; t++)
#pragma unroll
        for (int n = 0; n <= 1; n++) {
          const fp_type tmp = map[FLATTENED_3D(u0 + n, v0 + t, w0 + i, map_index_x, map_index_xy)];
          value += pu[n] * pv[t] * pw[i] * tmp;
        }
    return value;
  }
} // namespace mudock
