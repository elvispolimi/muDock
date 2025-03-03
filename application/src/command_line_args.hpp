#pragma once

#include <filesystem>
#include <mudock/knobs.hpp>
#include <vector>

static constexpr auto use_cpu_conf = std::string_view{"CPP:CPU:0"};

struct command_line_arguments {
  std::filesystem::path protein_path    = std::filesystem::path{"protein.pdb"};
  std::vector<std::string> device_confs = {std::string{use_cpu_conf}};
  mudock::knobs knobs;
};
command_line_arguments parse_command_line_arguments(const int argc, char *argv[]);
