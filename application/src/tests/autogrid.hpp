#pragma once

#include <filesystem>
#include <fstream>
#include <mudock/format/autogrid.hpp>
#include <mudock/log.hpp>
#include <mudock/type_alias.hpp>
#include <sstream>
#include <stdexcept>
#include <string>

inline std::string read_dpf(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("Cannot open DPF: " + path.string());
  std::ostringstream contents;
  contents << input.rdbuf();
  return contents.str();
}

inline std::filesystem::path resolve_dpf_path(const std::filesystem::path& dpf,
                                              const std::string& reference) {
  const std::filesystem::path path{reference};
  return path.is_absolute() ? path : dpf.parent_path() / path;
}

static inline mudock::autodock_grid load_autogrid_map_dpf(const std::filesystem::path& dpf_path) {
  return mudock::load_autogrid_maps(dpf_path);
}

static inline mudock::fp_type load_autodock_score(const std::filesystem::path& dpf_path) {
  std::istringstream input(read_dpf(dpf_path));
  std::string line;
  while (std::getline(input, line)) {
    const auto marker = line.find("AUTODOCK_SCORE");
    if (marker == std::string::npos) continue;
    std::istringstream value(line.substr(marker + std::string("AUTODOCK_SCORE").size()));
    mudock::fp_type score{};
    if (value >> score) return score;
  }
  throw std::runtime_error("DPF has no AUTODOCK_SCORE reference: " + dpf_path.string());
}

static inline mudock::fp_type load_autodock_error_score(const std::filesystem::path& dpf_path) {
  std::istringstream input(read_dpf(dpf_path));
  std::string line;
  while (std::getline(input, line)) {
    const auto marker = line.find("AUTODOCK_ERROR_CORRECTION");
    if (marker == std::string::npos) continue;
    std::istringstream value(line.substr(marker + std::string("AUTODOCK_ERROR_CORRECTION").size()));
    mudock::fp_type score{};
    if (value >> score) return score;
  }
  return mudock::fp_type{0};
}

static inline std::filesystem::path get_ligand_path(const std::filesystem::path& dpf_path) {
  std::istringstream input(read_dpf(dpf_path));
  std::string line;
  while (std::getline(input, line)) {
    std::istringstream fields(line);
    std::string directive, path;
    if (fields >> directive >> path; directive == "move") return resolve_dpf_path(dpf_path, path);
  }
  throw std::runtime_error("DPF has no move ligand directive: " + dpf_path.string());
}
