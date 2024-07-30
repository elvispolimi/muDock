#pragma once

#include <iostream>
#include <mudock/chem.hpp>
#include <mudock/grid/mdindex.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <stdexcept>
#include <vector>

namespace mudock {

  class grid_map {
    std::vector<fp_type> grid_values;
    index3D index;

  public:
    grid_map(const index3D& npts): index(npts) {
      grid_values.resize(npts.size_x() * npts.size_y() * npts.size_z());
    }
    ~grid_map() = default;
    grid_map(grid_map&& other): index(other.index) { grid_values.swap(other.grid_values); }
    grid_map(const grid_map& other)      = delete;
    grid_map& operator=(grid_map&&)      = default;
    grid_map& operator=(const grid_map&) = delete;

    inline fp_type& at(const size_t x, const size_t y, const size_t z) {
      return grid_values[index.to1D(x, y, z)];
    }
  };

  class grid_atom_map: public grid_map {
    autodock_ff atom_type;

  public:
    grid_atom_map(const autodock_ff& type, const index3D& npts): grid_map(npts), atom_type(type) {}
    ~grid_atom_map() = default;
    grid_atom_map(grid_atom_map&& other): grid_map(std::move(other)) { atom_type = other.atom_type; }
    grid_atom_map(const grid_atom_map& other)      = delete;
    grid_atom_map& operator=(grid_atom_map&&)      = default;
    grid_atom_map& operator=(const grid_atom_map&) = delete;

    [[nodiscard]] inline auto get_atom_type() const { return atom_type; }
    // Check this
    size_t is_hbonder{0};
  };

  std::vector<grid_atom_map> generate_atom_grid_maps(dynamic_molecule&);
  grid_map generate_electrostatic_grid_map(dynamic_molecule&);
  grid_map generate_desolvation_grid_map(dynamic_molecule&);
} // namespace mudock
