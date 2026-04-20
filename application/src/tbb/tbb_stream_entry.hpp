#pragma once

#include "../command_line_args.hpp"

#include <mudock/mpi_implementation/byte_range.hpp>

#include <optional>

int run_tbb_stream_entry(const command_line_arguments& args,
                         std::optional<mudock::byte_range> range = std::nullopt,
                         std::optional<int> rank = std::nullopt);
