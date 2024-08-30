#pragma once

#include <cmath>
#include <mudock/grid/mdvector.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  // This class model 3D space grid, where each pixel is a floating point. It can be accessed using the
  // coordinates rather than indexes
  class space_grid: public md_vector<fp_type, 3> {
    fp_type _inv_resolution = 2;

    // utility function that compute the flat index of the given point
    [[nodiscard]] std::size_t get_index(const point3D& p) const {
      return md_index<3>::to1D(static_cast<std::size_t>((p.x - _min.x) * _inv_resolution),
                               static_cast<std::size_t>((p.y - _min.y) * _inv_resolution),
                               static_cast<std::size_t>((p.z - _min.z) * _inv_resolution));
    }

  public:
    // information about the 3D space that we are representing
    point3D _min, _max, _center;

    // space grid constructors that try to figure out the space that we need to model
    inline space_grid(): md_vector<fp_type, 3>(1, 1, 1), _min(0, 0, 0), _max(0, 0, 0), _center(0, 0, 0) {}
    inline space_grid(const point3D min, const point3D max, const fp_type resolution)
        : md_vector<fp_type, 3>(std::ceil(std::abs(min.x - max.x) / resolution),
                                std::ceil(std::abs(min.y - max.y) / resolution),
                                std::ceil(std::abs(min.z - max.z) / resolution)),
          _inv_resolution(fp_type{1} / resolution),
          _min(min),
          _max(max),
          _center((max.x + min.x / fp_type{2}) + min.x,
                  (max.y + min.y / fp_type{2}) + min.y,
                  (max.z + min.z / fp_type{2}) + min.z) {}

    // function to check if the point fall inside the space grid
    inline bool is_outside(point3D p) {
      return (p.x < _min.x) || (p.x > _max.x) || (p.y < _min.y) || (p.y > _max.y) || (p.z < _min.z) ||
             (p.z > _max.z);
    }

    // function to access data of the space grid w/out checking if the point is actually inside
    inline fp_type& get(point3D p) { return md_vector<fp_type, 3>::_data[get_index(p)]; }
    inline const fp_type& get(point3D p) const { return md_vector<fp_type, 3>::_data[get_index(p)]; }
  };

} // namespace mudock
