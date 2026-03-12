#include <mpi.h>

#include <mudock/mpi_implementation/utilities.hpp>
#include <mudock/tbb_implementation/tbb_pipeline.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>

#include "../command_line_args.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <utility>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0, nranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  const auto args = parse_command_line_arguments(argc, argv);

  MUDOCK_MARKER_INIT;

  if (rank == 0) {
    mudock::info("Running with ", nranks, " MPI processes.");
  }

  {
    std::ifstream probe(args.ligand_path, std::ios::binary);
    if (!probe) {
      mudock::error("[rank ", rank, "] Cannot open ligand file: ", args.ligand_path);
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }

  const auto in_format = mudock::parse_supported_format(args.ligand_path);
  if (in_format != mudock::supported_format::ADTMOL2) {
    if (rank == 0) {
      mudock::error("MPI implementation currently supports only ADTMOL2 input format.");
    }
    MPI_Abort(MPI_COMM_WORLD, 3);
  }

  // Each rank reads the protein
  auto protein = std::make_shared<mudock::dynamic_molecule>(
      mudock::parser<mudock::dynamic_molecule>(args.protein_path));

  // Each rank computes the range of the input file it has to process    
  const auto range =
    mudock::mpi_splitter_bcast(args.ligand_path, rank, nranks);

  const size_t begin = range.first;
  const size_t end   = range.second;

  // print range per rank (in bytes)
  mudock::info("rank ", rank, ": [", begin, ", ", end, "]\n");

  //  compute our slab
  const bool did_work = (begin < end);
  if (did_work) {

    std::ifstream in(args.ligand_path, std::ios::binary);
    if (!in) {
      mudock::error("[rank ", rank, "] Can't open input file ", args.ligand_path);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    in.seekg(static_cast<std::streamoff>(begin), std::ios::beg);

    mudock::genetic_adt_pipeline pipe{protein};
    mudock::run_tbb_pipeline(in, args.device_confs, args.knobs, pipe, end);

  } else {
    mudock::info("[rank ", rank, "] No work (empty range).");
  }

  MUDOCK_MARKER_CLOSE;
  MPI_Finalize();

  // if we reach this statement we completed successfully the run
  mudock::info("All Done!");
  return EXIT_SUCCESS;
}
