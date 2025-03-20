#pragma once

#include <cassert>
#include <cstdlib>
#include <format>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct GridMap {
  int x_size, y_size, z_size;
  float origin[3]; // Origin of the grid (x, y, z)
  float spacing;   // Spacing between grid points (dx, dy, dz)
  std::vector<std::vector<std::vector<float>>> data;

  GridMap(int x, int y, int z)
      : x_size(x),
        y_size(y),
        z_size(z),
        data(x + 1, std::vector<std::vector<float>>(y + 1, std::vector<float>(z + 1, 0.0))) {}

  // Function to navigate the grid and get the value at (x, y, z)
  float get(const int x, const int y, const int z) const {
    if (x >= 0 && x < x_size && y >= 0 && y < y_size && z >= 0 && z < z_size) {
      return data[x][y][z];
    }
    throw std::runtime_error(std::format("Failed to access value at {} {} {}", x, y, z));
  }

  void set_map(const int x, const int y, const int z, const float val) {
    assert(x <= x_size && y <= y_size && z <= z_size);
    data[x][y][z] = val;
  }
};

// Function to read AutoDock map file and load it into GridMap
static inline GridMap loadGridMap(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + filename);
  }

  std::string line;
  int x = 0, y = 0, z = 0;
  float o_x = 0, o_y = 0, o_z = 0;
  float spacing = 0;

  // Read the header information (assuming standard AutoGrid format)
  while (std::getline(file, line)) {
    if (line.find("NELEMENTS") != std::string::npos) {
      // Read grid sizes and spacing
      sscanf(line.c_str(), "NELEMENTS %d %d %d", &x, &y, &z);
    } else if (line.find("SPACING") != std::string::npos) {
      sscanf(line.c_str(), "SPACING %f", &spacing);
    } else if (line.find("CENTER") != std::string::npos) {
      sscanf(line.c_str(), "CENTER %f %f %f", &o_x, &o_y, &o_z);
      break;
    }
  }

  GridMap gridMap{x, y, z};
  gridMap.origin[0] = o_x;
  gridMap.origin[1] = o_y;
  gridMap.origin[2] = o_z;
  gridMap.spacing   = spacing;

  // Read the grid data
  for (int v = 0; v <= z; ++v) {
    for (int t = 0; t <= y; ++t) {
      for (int u = 0; u <= x; u++) {
        std::getline(file, line);
        std::istringstream ss(line);
        float value;
        ss >> value;
        gridMap.set_map(u, t, v, value);
      }
    }
  }
  std::getline(file, line);
  assert(line.empty());

  file.close();
  return gridMap;
}
