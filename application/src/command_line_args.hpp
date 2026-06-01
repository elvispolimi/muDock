#pragma once

#include <filesystem>
#include <mudock/compute/algorithm.hpp>
#include <mudock/knobs.hpp>
#include <optional>
#include <string_view>
#include <vector>

namespace mudock {
    enum class scoring_mode {
        standard,
        precomputed,
        quant
    };

    inline scoring_mode parse_scoring_mode(const std::string& mode_str) {
        if (mode_str == "STANDARD" || mode_str == "standard") return scoring_mode::standard;
        if (mode_str == "PRECOMPUTED" || mode_str == "precomputed") return scoring_mode::precomputed;
        if (mode_str == "QUANT" || mode_str == "quant") return scoring_mode::quant;
        
        throw std::invalid_argument("Invalid scoring mode provided: " + mode_str);
    }
}

static constexpr auto use_cpu_conf = std::string_view{"CPP:CPU:0"};

struct command_line_arguments {
  std::filesystem::path protein_path    = std::filesystem::path{"protein.pdb"};
  std::filesystem::path ligand_path     = std::filesystem::path{"ligand.mol2"};
  std::vector<std::string> device_confs = {std::string{use_cpu_conf}};
  std::optional<double> time_limit_sec  = std::nullopt;
  std::optional<double> observer        = std::nullopt;
  mudock::search_algorithm search       = mudock::search_algorithm::NONE;
  mudock::scoring_function scoring      = mudock::scoring_function::ADT;
  std::string pipeline_mode = "STANDARD";
  mudock::scoring_mode score_mode       = mudock::scoring_mode::standard;
  bool score_only = false;
  mudock::knobs knobs;
};
command_line_arguments parse_command_line_arguments(const int argc, char *argv[]);
