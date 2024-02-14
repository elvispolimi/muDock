#pragma once

#include <filesystem>

struct command_line_arguments {
  std::filesystem::path protein_path = std::filesystem::path{"protein.pdb"};
};
command_line_arguments parse_command_line_arguments(const int argc, char *argv[]);
