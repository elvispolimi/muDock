#include "command_line_args.hpp"

#include <cstdlib>
#include <fstream>
#include <limits>
#include <memory>
#include <mudock/compute/pipeline_selector.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mpi_implementation/byte_range.hpp>
#include <mudock/mudock.hpp>
#include <mudock/tbb_implementation/tbb_pipeline.hpp>
#include <optional>

#ifdef MUDOCK_USE_MPI
  #include <mpi.h>
  #include <mudock/mpi_implementation/distributed_ranges.hpp>
#endif

namespace {
  int run_selected_pipeline(std::istream& in,
                            const mudock::supported_format format,
                            const command_line_arguments& args,
                            std::shared_ptr<mudock::dynamic_molecule> protein,
                            const std::uint64_t range_end) {
    mudock::info("Pipeline selection: search=", to_string(args.search), ", score=", to_string(args.scoring));

    dispatch_selected_pipeline(args.search,
                               args.scoring,
                               format,
                               [&]<typename pipeline_t>(const auto, const auto) {
                                 pipeline_t pipe{protein};
                                 constexpr_switch<0, mudock::get_num_supported_format(), 1>(
                                     [&](const auto format_index) {
                                       constexpr auto selected_format =
                                           static_cast<mudock::supported_format>(decltype(format_index)::value);
                                       mudock::run_tbb_pipeline<selected_format>(in,
                                                                                 args.device_confs,
                                                                                 args.knobs,
                                                                                 pipe,
                                                                                 range_end,
                                                                                 args.time_limit_sec,
                                                                                 args.observer);
                                     },
                                     format);
                               });

    return EXIT_SUCCESS;
  }
} // namespace

int main(int argc, char** argv) {
  const auto args                         = parse_command_line_arguments(argc, argv);
  std::optional<mudock::byte_range> range = std::nullopt;
  std::optional<int> rank                 = std::nullopt;
  const auto in_format                    = mudock::parse_supported_format(args.ligand_path);

#ifdef MUDOCK_USE_MPI
  MPI_Init(&argc, &argv);

  int mpi_rank = 0, nranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  rank = mpi_rank;

  if (rank == 0) {
    mudock::info("Running with ", nranks, " MPI processes.");
  }

  {
    std::ifstream probe(args.ligand_path, std::ios::binary);
    if (!probe) {
      mudock::error("[rank ", *rank, "] Cannot open ligand file: ", args.ligand_path);
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }

  constexpr_switch<0, mudock::get_num_supported_format(), 1>(
      [&](const auto format_index) {
        constexpr auto selected_format =
            static_cast<mudock::supported_format>(decltype(format_index)::value);
        range = mudock::distribute_aligned_ranges<selected_format>(args.ligand_path, *rank, nranks);
      },
      in_format);
  mudock::info("rank ", *rank, ": [", range->begin, ", ", range->end, "]\n");
#endif

  MUDOCK_MARKER_INIT;

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
#ifdef MUDOCK_USE_MPI
    MPI_Finalize();
#endif
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
#ifdef MUDOCK_USE_MPI
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
  }

  in.seekg(static_cast<std::streamoff>(effective_range.begin), std::ios::beg);

  const int status = run_selected_pipeline(in, in_format, args, std::move(protein), effective_range.end);
  MUDOCK_MARKER_CLOSE;
  if (status == EXIT_SUCCESS) {
    mudock::info("All Done!");
  }

#ifdef MUDOCK_USE_MPI
  MPI_Finalize();
#endif

  return status;
}
