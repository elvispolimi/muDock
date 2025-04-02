#pragma once

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <mudock/chem/grid_const.hpp>
#include <mudock/grid/grid_map.hpp>
#include <mudock/grid/mdindex.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <sstream>
#include <string>

struct autogrid_tokens {
  static constexpr auto ELEMENTS_TOKEN = "NELEMENTS";
  static constexpr auto SPACING_TOKEN  = "SPACING";
  static constexpr auto CENTER_TOKEN   = "CENTER";
};

// Function to read AutoDock map file and load it into GridMap
static inline mudock::grid_map load_autogrid_map(const std::string &filename) {
  const auto desc = read_from_stream(std::ifstream(filename));
  std::stringstream desc_s{desc};

  std::string line;
  mudock::index3D sizes;
  mudock::point3D center;

  // Read the header information (assuming standard AutoGrid format)
  while (std::getline(desc_s, line)) {
    std::istringstream stream(line);
    std::string _;
    if (line.find(autogrid_tokens::ELEMENTS_TOKEN) != std::string::npos) {
      int x, y, z;
      stream >> _ >> x >> y >> z;
      sizes = mudock::index3D{x + 1, y + 1, z + 1};
    } else if (line.find(autogrid_tokens::SPACING_TOKEN) != std::string::npos) {
      mudock::fp_type spacing;
      stream >> _ >> spacing;
      assert(spacing == mudock::grid_spacing);
    } else if (line.find(autogrid_tokens::CENTER_TOKEN) != std::string::npos) {
      mudock::fp_type x, y, z;
      stream >> _ >> x >> y >> z;
      center = mudock::point3D{x, y, z};
      // Assume that data are given in a certain order
      break;
    }
  }

  mudock::grid_map grid_map{sizes, center};

  // Read the grid data
  for (int v = 0; v < sizes.size_z(); ++v) {
    for (int t = 0; t < sizes.size_y(); ++t) {
      for (int u = 0; u < sizes.size_x(); ++u) {
        std::getline(desc_s, line);
        std::istringstream ss(line);
        mudock::fp_type value;
        ss >> value;
        grid_map.at(u, t, v) = value;
      }
    }
  }
  std::getline(desc_s, line);
  assert(line.empty());

  return grid_map;
}
