#include "../tbb/tbb_stream_entry.hpp"

#include <mpi.h>

#include <mudock/mpi_implementation/distributed_ranges.hpp>
#include <mudock/mudock.hpp>

#include "../command_line_args.hpp"

#include <fstream>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0, nranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  const auto args = parse_command_line_arguments(argc, argv);

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
  mudock::byte_range range{};
  constexpr_switch<0, mudock::get_num_supported_format(), 1>(
      [&](const auto format_index) {
        constexpr auto format = static_cast<mudock::supported_format>(format_index());
        if constexpr (format == mudock::supported_format::ADTMOL2) {
          range = mudock::distribute_aligned_ranges<format>(args.ligand_path, rank, nranks);
        }
      },
      in_format);

  // print range per rank (in bytes)
  mudock::info("rank ", rank, ": [", range.begin, ", ", range.end, "]\n");

  const int status = run_tbb_stream_entry(args, range, rank);
  MPI_Finalize();
  return status;
}
