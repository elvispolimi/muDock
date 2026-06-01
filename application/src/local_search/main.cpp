#include "command_line_args.hpp"

#include <cstdlib>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>

#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>
#include <mudock/mpi_implementation/byte_range.hpp>
#include <mudock/tbb_implementation/tbb_pipeline.hpp>

int main(int argc, char** argv) {
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                                                    //
  //  WARNING: to ensure a stan-alone application of local search, some knobs are forced in the local_search_pipeline,  //
  //           like population_number = 1 and num_generations = 1                                                       //
  //                                                                                                                    //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  const auto args = parse_command_line_arguments(argc, argv);
  std::optional<mudock::byte_range> range = std::nullopt;
  std::optional<int> rank                 = std::nullopt;
  const auto in_format                    = mudock::parse_supported_format(args.ligand_path);

  if (in_format != mudock::supported_format::ADTMOL2) {
    mudock::error("Stream implementation currently supports only ADTMOL2 input format.");
    return 1;
  }

  MUDOCK_MARKER_INIT;

  mudock::info("Reading and parsing protein ", args.protein_path, " ...");
  auto protein = std::make_shared<mudock::dynamic_molecule>(mudock::parser<mudock::dynamic_molecule>(args.protein_path));

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

  mudock::ls_adt_adadelta_pipeline pipe{protein};
          mudock::run_tbb_pipeline<mudock::supported_format::ADTMOL2>(
              in, args.device_confs, args.knobs, pipe, effective_range.end, args.time_limit_sec, args.observer);
  MUDOCK_MARKER_CLOSE;
  mudock::info("All Done!");

  return EXIT_SUCCESS;
}
