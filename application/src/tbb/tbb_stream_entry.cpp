#include "tbb_stream_entry.hpp"

#include <cstdlib>
#include <fstream>
#include <limits>
#include <memory>

#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>
#include <mudock/tbb_implementation/tbb_pipeline.hpp>

int run_tbb_stream_entry(const command_line_arguments& args,
                         std::optional<mudock::byte_range> range,
                         std::optional<int> rank) {
  MUDOCK_MARKER_INIT;

  const auto in_format = mudock::parse_supported_format(args.ligand_path);
  if (in_format != mudock::supported_format::ADTMOL2 && in_format != mudock::supported_format::MOL2) {
    if (rank) {
      mudock::error("MPI implementation currently supports only ADTMOL2 and MOL2 input formats.");
    } else {
      mudock::error("TBB implementation currently supports only ADTMOL2 and MOL2 input formats.");
    }
    return 1;
  }

  mudock::info("Reading and parsing protein ", args.protein_path, " ...");
  auto protein =
      std::make_shared<mudock::dynamic_molecule>(mudock::parser<mudock::dynamic_molecule>(args.protein_path));

  mudock::info("Reading ligand ", args.ligand_path, " ...");
  std::ifstream in(args.ligand_path, std::ios::binary);
  if (!in) {
    if (rank) {
      mudock::error("[rank ", *rank, "] Can't open input file ", args.ligand_path);
    } else {
      mudock::error("Can't open input file ", args.ligand_path);
    }
    return 1;
  }

  const auto effective_range =
      range.value_or(mudock::byte_range{0, std::numeric_limits<std::uint64_t>::max()});
  if (effective_range.empty()) {
    if (rank) {
      mudock::info("[rank ", *rank, "] No work (empty range).");
    }
    MUDOCK_MARKER_CLOSE;
    mudock::info("All Done!");
    return EXIT_SUCCESS;
  }

  in.seekg(static_cast<std::streamoff>(effective_range.begin), std::ios::beg);

  mudock::genetic_adt_pipeline pipe{protein};
  constexpr_switch<0, mudock::get_num_supported_format(), 1>(
      [&](const auto format_index) {
        constexpr auto format = static_cast<mudock::supported_format>(format_index());
        if constexpr (format == mudock::supported_format::ADTMOL2 || format == mudock::supported_format::MOL2) {
          mudock::run_tbb_pipeline<format>(
              in, args.device_confs, args.knobs, pipe, effective_range.end, args.time_limit_sec, args.observer);
        }
      },
      in_format);

  MUDOCK_MARKER_CLOSE;
  mudock::info("All Done!");
  return EXIT_SUCCESS;
}
