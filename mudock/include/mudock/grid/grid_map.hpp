#pragma once

#include <cmath>
#include <cstddef>
#include <iostream>
#include <mudock/chem.hpp>
#include <mudock/grid/mdindex.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <stdexcept>
#include <vector>

namespace mudock {
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
      epsilon = 1.0L;
    }
    return epsilon;
  }

  template<class T, class index_type>
    requires is_index<index_type>
  class grid {
    std::vector<T> grid_values;

  public:
    grid(const index_type _index): grid_values(_index.get_dim()), index(_index) {};
    ~grid()                       = default;
    grid(grid&& other)            = default;
    grid(const grid& other)       = delete;
    grid& operator=(grid&& other) = default;
    grid& operator=(const grid&)  = delete;

    const index_type index;

    template<typename... Indexes>
    [[nodiscard]] inline T& at(Indexes... indexes) {
      return grid_values[index.to1D(indexes...)];
    }

    template<typename... Indexes>
    [[nodiscard]] inline T at(Indexes... indexes) const {
      return grid_values[index.to1D(indexes...)];
    }
  };

  class grid_map: public grid<fp_type, index3D> {
  public:
    const point3D minimum, maximum, center;
    grid_map(const index3D& npts, const point3D& min, const point3D& max)
        : grid(npts),
          minimum(min),
          maximum(max),
          center{(maximum.x - minimum.x) / 2 + minimum.x,
                 (maximum.y - minimum.y) / 2 + minimum.y,
                 (maximum.z - minimum.z) / 2 + minimum.z} {}
    ~grid_map()                           = default;
    grid_map(grid_map&& other)            = default;
    grid_map(const grid_map& other)       = delete;
    grid_map& operator=(grid_map&& other) = default;
    grid_map& operator=(const grid_map&)  = delete;

    [[nodiscard]] inline int outside_grid(const point3D& p) const {
      return p.x < minimum.x || p.x > maximum.x || p.y < minimum.y || p.y > maximum.y || p.z < minimum.z ||
             p.z > maximum.z;
    }

    [[nodiscard]] inline point3D get_index_from_coordinates(const point3D coord) const {
      return {(coord.x - minimum.x) / grid_spacing,
              (coord.y - minimum.y) / grid_spacing,
              (coord.z - minimum.z) / grid_spacing};
    }

    [[nodiscard]] inline fp_type& at(const point3D& p) { return grid::at(p.x, p.y, p.z); }

    [[nodiscard]] inline fp_type at(const point3D& p) const { return grid::at(p.x, p.y, p.z); }

    [[nodiscard]] inline fp_type& at(const std::size_t x, const std::size_t y, const std::size_t z) {
      return grid::at(x, y, z);
    }

    [[nodiscard]] inline fp_type at(const std::size_t x, const std::size_t y, const std::size_t z) const {
      return grid::at(x, y, z);
    }
  };

  class grid_atom_map: public grid_map {
    autodock_ff atom_type;

  public:
    grid_atom_map(const autodock_ff& type, const index3D& npts, const point3D& min, const point3D& max)
        : grid_map(npts, min, max), atom_type(type) {}
    ~grid_atom_map()                               = default;
    grid_atom_map(grid_atom_map&& other)           = default;
    grid_atom_map(const grid_atom_map& other)      = delete;
    grid_atom_map& operator=(grid_atom_map&&)      = default;
    grid_atom_map& operator=(const grid_atom_map&) = delete;

    [[nodiscard]] inline auto get_atom_type() const { return atom_type; }
    // Check this
    size_t is_hbonder{0};
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
  };

  [[nodiscard]] grid_atom_mapper generate_atom_grid_maps(dynamic_molecule&);
  [[nodiscard]] grid_map generate_electrostatic_grid_map(dynamic_molecule&);
  [[nodiscard]] grid_map generate_desolvation_grid_map(dynamic_molecule&);
} // namespace mudock
