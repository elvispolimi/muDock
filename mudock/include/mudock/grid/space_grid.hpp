#pragma once

#include <cmath>
#include <mudock/grid/mdindex.hpp>
#include <mudock/grid/mdspan.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/type_alias.hpp>
#include <type_traits>

namespace mudock {

  // This class model 3D space grid, where each pixel is a floating point. It can be accessed using the
  // coordinates rather than indexes
  // FIXME N==3 for the time being, requires making point<fp_type, n>.hpp operations parametric on the size of the point
  template<typename C, typename T, std::size_t n>
  // requires((std::is_same<C, md_span<T, n>>::value || std::is_same<C, md_container<T, n>>::value) && n == 3)
  class space_grid_t {
    C md_data;

    // utility function that compute the flat index of the given point
    [[nodiscard]] std::size_t get_index(const point<fp_type, n>& p) const {
      // return md_index<n>::to1D(static_cast<std::size_t>((p.x - _min.x) * _inv_resolution),
      //                          static_cast<std::size_t>((p.y - _min.y) * _inv_resolution),
      //                          static_cast<std::size_t>((p.z - _min.z) * _inv_resolution));
      return md_index<n>::to1D((p - _min).truncate() * _inv_resolution);
    }

  public:
    // information about the 3D space that we are representing
    fp_type _inv_resolution = 2;
    point<fp_type, n> _min, _max, _center;

    inline space_grid_t(): md_data() {}
    space_grid_t(const point<fp_type, n> min,
                 const point<fp_type, n> max,
                 const point<fp_type, n> center,
                 const fp_type resolution,
                 const C data)
        // requires(std::is_same<C, md_span<T, n>>::value)
        : md_data(data), _inv_resolution(fp_type{1} / resolution), _min(min), _max(max), _center(center) {}
    template<class... Y>
    space_grid_t(const point<fp_type, n> min,
                 const point<fp_type, n> max,
                 const point<fp_type, n> center,
                 const fp_type resolution,
                 Y&&... sizes)
        // requires(std::is_same<C, md_container<T, n>>::value && sizeof...(sizes) == n)
        : md_data(sizes...),
          _inv_resolution(fp_type{1} / resolution),
          _min(min),
          _max(max),
          _center(center) {}

    // function to check if the point fall inside the space grid
    [[nodiscard]] inline bool is_outside(point<fp_type, n> p) {
      // return (p.x < _min.x) || (p.x > _max.x) || (p.y < _min.y) || (p.y > _max.y) || (p.z < _min.z) ||
      //        (p.z > _max.z);
      return (p < _min) || (p > _max);
    }

    template<std::size_t index>
    [[nodiscard]] std::size_t size() const {
      return md_data.template size<index>();
    }

    // function to convert the index to the coordinate of the point
    template<typename... I>
    [[nodiscard]] inline auto to_coord(I... sizes) const
      requires(sizeof...(sizes) == n)
    {
      return (point<T, n>{sizes...} * _inv_resolution) + _min;
    }

    [[nodiscard]] T x() const { return md_data.size_x(); }
    [[nodiscard]] T y() const
      requires(n > 1)
    {
      return md_data.size_y();
    }
    [[nodiscard]] T z() const
      requires(n > 2)
    {
      return md_data.size_z();
    }

    // function to access data of the space grid w/out checking if the point is actually inside
    [[nodiscard]] inline fp_type& get(point<fp_type, n> p) { return md_data._data[get_index(p)]; }
    [[nodiscard]] inline const fp_type& get(point<fp_type, n> p) const { return md_data._data[get_index(p)]; }

    // forward functions to access data using indexes instead of coordinates
    template<class... Y>
    [[nodiscard]] inline fp_type& get(Y&&... indexes) {
      return md_data.get(indexes...);
    }
    template<class... Y>
    [[nodiscard]] inline const fp_type& get(Y&&... indexes) const {
      return md_data.get(indexes...);
    }

    [[nodiscard]] inline const md_index<n>& get_space_index() const { return md_data; }
    [[nodiscard]] inline T* data() const { return md_data.data(); };
  };

  template<typename T>
  using space_grid_view = space_grid_t<md_span<T, 3>, T, 3>;
  // using space_grid_view_const = space_grid_t<md_span<const fp_type, 3>, const fp_type, 3>;
  using space_grid = space_grid_t<md_container<std::vector<fp_type>, 3>, std::vector<fp_type>, 3>;

} // namespace mudock
