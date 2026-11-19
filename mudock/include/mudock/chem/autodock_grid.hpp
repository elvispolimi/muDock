#pragma once

#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/grid/mdindex.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/grid/space_grid.hpp>
#include <mudock/type_alias.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace mudock {
  struct autodock_grid {
    autodock_grid(const md_index<3> dimensions, const point3D min, const point3D max)
        : index(dimensions),
          data(index.size_x(), index.size_y(), index.size_z(), num_autodock_grids()),
          _resolution{(max.x() - min.x()) /
                      static_cast<fp_type>(index.size_x() > 1 ? index.size_x() - 1 : 1)},
          _min(min),
          _max(max),
          _center{max.difference(min).divide({fp_type{2}}).add(min)},
          present{} {
      if (index.size_x() < 2 || index.size_y() < 2 || index.size_z() < 2)
        throw std::invalid_argument("AutoGrid dimensions must be at least 2 in every axis");
      const auto sy = (max.y() - min.y()) / static_cast<fp_type>(index.size_y() - 1);
      const auto sz = (max.z() - min.z()) / static_cast<fp_type>(index.size_z() - 1);
      const auto relative_tolerance = static_cast<fp_type>(0.00001);
      const auto tolerance = std::max(relative_tolerance, std::abs(_resolution) * relative_tolerance);
      if (!std::isfinite(_resolution) || !std::isfinite(sy) || !std::isfinite(sz) || _resolution <= 0 || sy <= 0 ||
          sz <= 0 || std::abs(sy - _resolution) > tolerance || std::abs(sz - _resolution) > tolerance)
        throw std::invalid_argument("AutoGrid spacing must be equal on all axes");
    }

    autodock_grid() = default;

    inline auto get_atom_map(const autodock_grid_type type) {
      return space_grid_view<fp_type>{_min,
                                      _max,
                                      _center,
                                      _resolution,
                                      data.get_slice(md_index<4>{index.size_x(),
                                                                 index.size_y(),
                                                                 index.size_z(),
                                                                 static_cast<int>(type)},
                                                     index)};
    }
    inline auto get_atom_map(const autodock_grid_type type) const {
      return space_grid_view<const fp_type>{_min,
                                            _max,
                                            _center,
                                            _resolution,
                                            data.get_slice(md_index<4>{index.size_x(),
                                                                       index.size_y(),
                                                                       index.size_z(),
                                                                       static_cast<int>(type)},
                                                           index)};
    }

    [[nodiscard]] inline auto get_map_flat_size() const { return index.flat_size(); }
    [[nodiscard]] inline const auto* get_maps_pointer() const { return data.data(); }
    [[nodiscard]] inline auto* get_min_p() const { return _min.data(); }
    [[nodiscard]] inline auto* get_max_p() const { return _max.data(); }
    [[nodiscard]] inline auto* get_center_p() const { return _center.data(); }
    [[nodiscard]] inline auto get_center() const { return _center; }
    [[nodiscard]] inline auto get_min() const { return _min; }
    [[nodiscard]] inline auto get_max() const { return _max; }
    [[nodiscard]] inline auto get_resolution() const { return _resolution; }
    [[nodiscard]] inline auto get_size_x() const { return index.size_x(); }
    [[nodiscard]] inline auto get_size_xy() const { return index.size_xy(); }
    [[nodiscard]] inline auto get_size_xyz() const { return index.flat_size(); }
    [[nodiscard]] inline bool has_map(const autodock_grid_type type) const {
      return present[static_cast<std::size_t>(type)];
    }
    inline void mark_map_present(const autodock_grid_type type) {
      present[static_cast<std::size_t>(type)] = true;
    }

    md_index<3> index;
    md_container<std::vector<fp_type>, 4> data;

  private:
    fp_type _resolution = fp_type{0.5};
    point3D _min, _max, _center;
    std::array<bool, num_autodock_grids()> present{};
  };
} // namespace mudock
