#pragma once

#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/grid/mdindex.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/grid/space_grid.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
using namespace std::placeholders;

struct fld_tokens {
  static constexpr auto FILE_REGEX     = "(file=([^\\s]+))";
  static constexpr auto FILE_TOKEN     = "file=";
  static constexpr auto VARIABLE_TOKEN = "variable";
  static constexpr auto LABEL_TOKEN    = "label";
  static constexpr auto ELETRO_TOKEN   = "Electrostatics";
  static constexpr auto DESOLV_TOKEN   = "Desolvation";
  static constexpr auto XYZ_TOKEN      = "xyz";
  static constexpr auto X_TOKEN        = "dim1";
  static constexpr auto Y_TOKEN        = "dim2";
  static constexpr auto Z_TOKEN        = "dim3";
};

struct dpf_tokens {
  static constexpr auto MAP_TOKEN    = "map ";
  static constexpr auto ELEC_TOKEN   = "elecmap";
  static constexpr auto DESOLV_TOKEN = "dsolvmap";
  static constexpr auto MOVE_TOKEN   = "move";
  static constexpr auto FLD_TOKEN    = "fld";
};

static constexpr auto SCORE_TOKEN       = "AUTODOCK_SCORE";
static constexpr auto ERROR_SCORE_TOKEN = "AUTODOCK_ERROR_CORRECTION";

static inline std::filesystem::path resolve_autodock_reference_path(const std::filesystem::path& owner_path,
                                                                    const std::string& raw_path) {
  const std::filesystem::path parsed_path{raw_path};
  if (std::filesystem::exists(parsed_path))
    return parsed_path;

  const auto fallback = owner_path.parent_path() / parsed_path.filename();
  if (std::filesystem::exists(fallback))
    return fallback;

  return parsed_path;
}

static inline mudock::autodock_grid load_autogrid_map_fld(const std::string& fld_path) {
  mudock::info("Reading and parsing FLD file ", fld_path, " for autogrid maps ...");

  std::string line;
  mudock::md_index<3> sizes;
  mudock::point3D min, max;
  std::array<std::filesystem::path, mudock::num_autodock_grids()> grids_filepath;
  int label{0};
  // FLD labels include affinity maps plus electrostatics and desolvation.
  std::array<mudock::autodock_grid_type, mudock::num_autodock_grids()> variables;
  const auto fld_desc = read_from_stream(std::ifstream(fld_path));
  const std::filesystem::path fld_fs_path{fld_path};
  std::stringstream fld_desc_s{fld_desc};

  int x{0}, y{0}, z{0};
  while (std::getline(fld_desc_s, line)) {
    std::istringstream stream(line);
    if (line.find(fld_tokens::X_TOKEN) != std::string::npos) {
      size_t equal_pos = line.find('=');
      if (equal_pos != std::string::npos) {
        x = std::stoi(line.substr(equal_pos + 1));
      }
    } else if (line.find(fld_tokens::Y_TOKEN) != std::string::npos) {
      size_t equal_pos = line.find('=');
      if (equal_pos != std::string::npos) {
        y = std::stoi(line.substr(equal_pos + 1));
      }
    } else if (line.find(fld_tokens::Z_TOKEN) != std::string::npos) {
      size_t equal_pos = line.find('=');
      if (equal_pos != std::string::npos) {
        z = std::stoi(line.substr(equal_pos + 1));
      }
    } else if (line.find(fld_tokens::XYZ_TOKEN) != std::string::npos) {
      std::smatch match;
      if (std::regex_search(line, match, std::regex(fld_tokens::FILE_REGEX))) {
        const auto xyz_path = resolve_autodock_reference_path(fld_fs_path, match[2]);
        const auto xyz_desc = read_from_stream(std::ifstream(xyz_path));

        std::stringstream xyz_desc_s{xyz_desc};

        mudock::fp_type min_x, max_x, min_y, max_y, min_z, max_z;
        std::getline(xyz_desc_s, line);
        std::stringstream xyz_stream{line};
        xyz_stream >> min_x >> max_x;

        std::getline(xyz_desc_s, line);
        xyz_stream.clear();
        xyz_stream.str(line);
        xyz_stream >> min_y >> max_y;

        std::getline(xyz_desc_s, line);
        xyz_stream.clear();
        xyz_stream.str(line);
        xyz_stream >> min_z >> max_z;

        std::getline(xyz_desc_s, line);
        assert(line.empty() && "XYZ file lenght is wrong");

        min = {min_x, min_y, min_z};
        max = {max_x, max_y, max_z};
      }
    } else if (line.find(fld_tokens::LABEL_TOKEN) != std::string::npos) {
      size_t equal_pos = line.find('=');
      size_t dash_pos  = line.find('-');

      if (equal_pos != std::string::npos && dash_pos != std::string::npos && equal_pos < dash_pos) {
        const std::string result = line.substr(equal_pos + 1, dash_pos - equal_pos - 1);
        variables[label]         = mudock::parse_map_symbol(result);
      }
      if (equal_pos != std::string::npos) {
        if (line.find(fld_tokens::ELETRO_TOKEN) != std::string::npos)
          variables[label] = mudock::autodock_grid_type::ELEC;
        else if (line.find(fld_tokens::DESOLV_TOKEN) != std::string::npos)
          variables[label] = mudock::autodock_grid_type::DESOLV;
      } else {
        throw std::runtime_error("Invalid label in fld");
      }
      ++label;
    } else if (line.find(fld_tokens::VARIABLE_TOKEN) != std::string::npos &&
               line.find(fld_tokens::FILE_TOKEN) != std::string::npos) {
      int id;

      std::istringstream v_stream(line);
      std::string token, map_path;
      v_stream >> token >> id >> map_path >> token;

      map_path = map_path.substr(std::strlen(fld_tokens::FILE_TOKEN));
      grids_filepath[static_cast<int>(variables[id - 1])] =
          resolve_autodock_reference_path(fld_fs_path, map_path);
    }
  }
  sizes = {x, y, z};

  mudock::autodock_grid adt_grid{min, max, mudock::grid_spacing};
  assert(adt_grid.index == static_cast<mudock::md_index<3>>(sizes));

  for (int map_index = 0; map_index < mudock::num_autodock_grids(); ++map_index) {
    if (grids_filepath[map_index].empty())
      continue;
    auto map = adt_grid.get_atom_map(static_cast<mudock::autodock_grid_type>(map_index));

    const auto map_desc = read_from_stream(std::ifstream(grids_filepath[map_index]));
    std::stringstream map_desc_s{map_desc};

    for (int index = 0; index < 6; ++index) std::getline(map_desc_s, line);

    // Read the grid data
    for (size_t v = 0; v < sizes.size_z(); ++v) {
      for (size_t t = 0; t < sizes.size_y(); ++t) {
        for (size_t u = 0; u < sizes.size_x(); ++u) {
          std::getline(map_desc_s, line);
          std::istringstream ss(line);
          mudock::fp_type value{0};
          ss >> value;
          map.get(u, t, v) = value;
        }
      }
    }
    std::getline(map_desc_s, line);
    assert(line.empty());
  }

  return adt_grid;
}

// Function to read AutoDock map file and load it into GridMap
static inline mudock::autodock_grid load_autogrid_map_dpf(const std::string& dpf_path) {
  mudock::info("Reading and parsing DPF file ", dpf_path, " for autogrid maps ...");
  const auto desc = read_from_stream(std::ifstream(dpf_path));
  const std::filesystem::path dpf_fs_path{dpf_path};
  std::stringstream desc_s{desc};

  std::string line;
  std::string fld_path, _;
  while (std::getline(desc_s, line)) {
    // Skip empty lines
    if (line.empty())
      continue;

    if (line.find(dpf_tokens::FLD_TOKEN) != std::string::npos) {
      std::stringstream ss{line};
      ss >> _ >> fld_path;
    }
  }

  return load_autogrid_map_fld(resolve_autodock_reference_path(dpf_fs_path, fld_path).string());
}

static inline mudock::fp_type load_autodock_score(const std::string& dpf_path) {
  mudock::info("Reading and parsing DPF file ", dpf_path, " for autodock error ...");
  const auto desc = read_from_stream(std::ifstream(dpf_path));
  std::stringstream desc_s{desc};

  std::string line;
  mudock::fp_type adt_score{0};
  // TODO check if you can get rid of this and use mudock autogrid maps
  while (std::getline(desc_s, line)) {
    // Skip empty lines
    if (line.empty())
      continue;

    if (line.find(SCORE_TOKEN) != std::string::npos) {
      std::stringstream ss{line};
      std::string _;
      ss >> _ >> _ >> adt_score;
      break;
    }
  }
  return adt_score;
}

static inline mudock::fp_type load_autodock_error_score(const std::string& dpf_path) {
  mudock::info("Reading and parsing DPF file ", dpf_path, " for autodock error score ...");
  const auto desc = read_from_stream(std::ifstream(dpf_path));
  std::stringstream desc_s{desc};

  std::string line;
  mudock::fp_type adt_error_score{0};
  // TODO check if you can get rid of this and use mudock autogrid maps
  while (std::getline(desc_s, line)) {
    // Skip empty lines
    if (line.empty())
      continue;

    if (line.find(ERROR_SCORE_TOKEN) != std::string::npos) {
      std::stringstream ss{line};
      std::string _;
      ss >> _ >> _ >> adt_error_score;
      break;
    }
  }
  return adt_error_score;
}

static inline std::string get_ligand_path(const std::string& dpf_path) {
  mudock::info("Reading and parsing DPF file ", dpf_path, " for ligand ...");
  const auto desc = read_from_stream(std::ifstream(dpf_path));
  const std::filesystem::path dpf_fs_path{dpf_path};
  std::stringstream desc_s{desc};

  std::string line;
  // TODO check if you can get rid of this and use mudock autogrid maps
  while (std::getline(desc_s, line)) {
    // Skip empty lines
    if (line.empty())
      continue;

    if (line.find(dpf_tokens::MOVE_TOKEN) != std::string::npos) {
      std::stringstream ss{line};
      std::string ligand_path, _;
      ss >> _ >> ligand_path;

      return resolve_autodock_reference_path(dpf_fs_path, ligand_path).string();
    }
  }
  throw std::runtime_error("Cannot file ligand path inside given PDF file");
}
