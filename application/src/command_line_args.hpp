#pragma once

#include "pipeline_selection.hpp"

#include <filesystem>
#include <mudock/knobs.hpp>
#include <optional>
#include <string_view>
#include <vector>

static constexpr auto use_cpu_conf = std::string_view{"CPP:CPU:0"};

struct command_line_arguments {
  std::filesystem::path protein_path    = std::filesystem::path{"protein.pdb"};
  std::filesystem::path ligand_path     = std::filesystem::path{"ligand.mol2"};
  std::vector<std::string> device_confs = {std::string{use_cpu_conf}};
  std::optional<double> time_limit_sec  = std::nullopt;
  std::optional<double> observer        = std::nullopt;
  search_algorithm search               = search_algorithm::NONE;
  scoring_function scoring             = scoring_function::ADT;
  mudock::knobs knobs;
};
command_line_arguments parse_command_line_arguments(const int argc, char *argv[]);
