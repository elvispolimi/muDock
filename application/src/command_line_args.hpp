#pragma once

#include <filesystem>
#include <mudock/compute/algorithm.hpp>
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
  mudock::search_algorithm search       = mudock::search_algorithm::NONE;
  mudock::scoring_function scoring      = mudock::scoring_function::ADT;
  mudock::knobs knobs;
};
command_line_arguments parse_command_line_arguments(const int argc, char *argv[]);
