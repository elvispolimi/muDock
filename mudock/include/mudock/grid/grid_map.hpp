#pragma once

#include "mudock/type_alias.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <mudock/chem.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/grid/mdindex.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <stdexcept>
#include <vector>

namespace mudock {
  // Utility for grid map creation

  inline point3D get_grid_minimum(const dynamic_molecule& receptor) {
    const fp_type receptor_min_x = std::ranges::min(receptor.get_x());
    const fp_type receptor_min_y = std::ranges::min(receptor.get_y());
    const fp_type receptor_min_z = std::ranges::min(receptor.get_z());
    return {std::floor((receptor_min_x - cutoff_distance - grid_spacing) * inv_spacing) / inv_spacing,
            std::floor((receptor_min_y - cutoff_distance - grid_spacing) * inv_spacing) / inv_spacing,
            std::floor((receptor_min_z - cutoff_distance - grid_spacing) * inv_spacing) / inv_spacing};
  }

  inline point3D get_grid_maximum(const dynamic_molecule& receptor) {
    const fp_type receptor_max_x = std::ranges::max(receptor.get_x());
    const fp_type receptor_max_y = std::ranges::max(receptor.get_y());
    const fp_type receptor_max_z = std::ranges::max(receptor.get_z());
    return {std::ceil((receptor_max_x + cutoff_distance + grid_spacing) * inv_spacing) / inv_spacing,
            std::ceil((receptor_max_y + cutoff_distance + grid_spacing) * inv_spacing) / inv_spacing,
            std::ceil((receptor_max_z + cutoff_distance + grid_spacing) * inv_spacing) / inv_spacing};
  }

  inline index3D get_grid_points_dim(const point3D& grid_minimum, const point3D& grid_maximum) {
    return {static_cast<int>((grid_maximum.x - grid_minimum.x) / grid_spacing),
            static_cast<int>((grid_maximum.y - grid_minimum.y) / grid_spacing),
            static_cast<int>((grid_maximum.z - grid_minimum.z) / grid_spacing)};
  }

  // TODO move it away from here
  inline fp_type calc_ddd_Mehler_Solmajer(fp_type distance) {
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
      epsilon = fp_type{1.0};
    }
    return epsilon;
  }

  template<class T, class index_type>
    requires is_index<index_type>
  class grid {
    std::vector<T> grid_values;

  public:
    grid(const index_type _index): grid_values(_index.get_dim()), index(_index){};
    ~grid()                       = default;
    grid(grid&& other)            = default;
    grid(const grid& other)       = default;
    grid& operator=(grid&& other) = default;
    grid& operator=(const grid&)  = default;

    const index_type index;

    [[nodiscard]] inline auto* data() const { return grid_values.data(); }

    template<typename... Indexes>
    [[nodiscard]] inline T& at(Indexes... indexes) {
      return grid_values[index.to1D(indexes...)];
    }

    template<typename... Indexes>
    [[nodiscard]] inline T at(Indexes... indexes) const {
      return grid_values[index.to1D(indexes...)];
    }

    void reset() { std::fill(grid_values.begin(), grid_values.end(), 0); }
  };

  class grid_map: public grid<fp_type, index3D> {
  public:
    const point3D minimum, maximum, center;
    const point3D minimum_coord, maximum_coord;

    grid_map(const dynamic_molecule& receptor)
        : grid(get_grid_points_dim(get_grid_minimum(receptor), get_grid_maximum(receptor))),
          minimum(get_grid_minimum(receptor)),
          maximum(get_grid_maximum(receptor)),
          center{(maximum.x - minimum.x) / 2 + minimum.x,
                 (maximum.y - minimum.y) / 2 + minimum.y,
                 (maximum.z - minimum.z) / 2 + minimum.z},
          minimum_coord({minimum.x + grid_spacing, minimum.y + grid_spacing, minimum.z + grid_spacing}),
          maximum_coord({maximum.x - grid_spacing, maximum.y - grid_spacing, maximum.z - grid_spacing}) {}
    ~grid_map()                           = default;
    grid_map(grid_map&& other)            = default;
    grid_map(const grid_map& other)       = default;
    grid_map& operator=(grid_map&& other) = delete;
    grid_map& operator=(const grid_map&)  = delete;

    [[nodiscard]] inline int outside_grid(const point3D& p) const {
      // TODO fix me to enable interpolation
      // Could be optimized
      return p.x < minimum_coord.x || p.x > maximum_coord.x || p.y < minimum_coord.y ||
             p.y > maximum_coord.y || p.z < minimum_coord.z || p.z > maximum_coord.z;
    }

    [[nodiscard]] inline point3D get_index_from_coordinates(const point3D coord) const {
      return {(coord.x - minimum.x) / grid_spacing,
              (coord.y - minimum.y) / grid_spacing,
              (coord.z - minimum.z) / grid_spacing};
    }

    [[nodiscard]] inline fp_type& at(const point3D& p) { return grid::at(p.x, p.y, p.z); }

    [[nodiscard]] inline fp_type at(const point3D& p) const { return grid::at(p.x, p.y, p.z); }

    [[nodiscard]] inline fp_type& at(const int x, const int y, const int z) { return grid::at(x, y, z); }

    [[nodiscard]] inline fp_type at(const int x, const int y, const int z) const { return grid::at(x, y, z); }

    [[nodiscard]] inline auto* data() const { return grid::data(); }
  };

  class grid_atom_map: public grid_map {
    autodock_ff atom_type;

  public:
    grid_atom_map(const autodock_ff& type, const dynamic_molecule& receptor)
        : grid_map(receptor), atom_type(type) {}
    ~grid_atom_map()                               = default;
    grid_atom_map(grid_atom_map&& other)           = default;
    grid_atom_map(const grid_atom_map& other)      = default;
    grid_atom_map& operator=(grid_atom_map&&)      = delete;
    grid_atom_map& operator=(const grid_atom_map&) = delete;

    [[nodiscard]] inline auto get_atom_type() const { return atom_type; }
    // Check this
    int is_hbonder{0};
  };

  // TODO check if we should put it together with also other maps
  class grid_atom_mapper {
    // TODO @Davide
    std::unordered_map<autodock_ff, grid_atom_map> grid_maps;

  public:
    grid_atom_mapper(std::vector<grid_atom_map>& maps) {
      for (auto& map: maps) { grid_maps.emplace(map.get_atom_type(), std::move(map)); }
    }

    [[nodiscard]] inline const grid_atom_map& get_atom_map(const autodock_ff& type) const {
      return grid_maps.at(type);
    }

    [[nodiscard]] inline const point3D& get_minimum() const { return grid_maps.begin()->second.minimum; }

    [[nodiscard]] inline const point3D& get_maximum() const { return grid_maps.begin()->second.maximum; }

    [[nodiscard]] inline const point3D& get_center() const { return grid_maps.begin()->second.center; }

    [[nodiscard]] inline const index3D& get_index() const { return grid_maps.begin()->second.index; }
  };

  [[nodiscard]] grid_atom_mapper generate_atom_grid_maps(dynamic_molecule&);
  [[nodiscard]] grid_map generate_electrostatic_grid_map(dynamic_molecule&);
  [[nodiscard]] grid_map generate_desolvation_grid_map(dynamic_molecule&);
} // namespace mudock
