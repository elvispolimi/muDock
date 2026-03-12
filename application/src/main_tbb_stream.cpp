#include "command_line_args.hpp"

#include <fstream>
#include <cstdlib>

#include <mudock/tbb_implementation/tbb_pipeline.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>

int main(int argc, char** argv) {
  const auto args = parse_command_line_arguments(argc, argv);

  MUDOCK_MARKER_INIT;

  const auto in_format = mudock::parse_supported_format(args.ligand_path);
  if (in_format != mudock::supported_format::ADTMOL2) {
    mudock::error("TBB implementation currently supports only ADTMOL2 input format.");
    return 1;
  }

  // read and parse the target protein
  mudock::info("Reading and parsing protein ", args.protein_path, " ...");
  auto protein =
      std::make_shared<mudock::dynamic_molecule>(mudock::parser<mudock::dynamic_molecule>(args.protein_path));

  mudock::info("Reading ligand ", args.ligand_path, " ...");
  std::ifstream in(args.ligand_path, std::ios::binary);
  if (!in) {
      mudock::error("Can't open input file ", args.ligand_path);
      return 1;
  }

  mudock::genetic_adt_pipeline pipe{protein};

  mudock::run_tbb_pipeline(in, args.device_confs, args.knobs, pipe);

  MUDOCK_MARKER_CLOSE;

  // if we reach this statement we completed successfully the run
  mudock::info("All Done!");
  return EXIT_SUCCESS;
}
