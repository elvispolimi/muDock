#include "command_line_args.hpp"

#include <boost/program_options.hpp>
#include <cstddef>
#include <iostream>
#include <optional>

command_line_arguments parse_command_line_arguments(const int argc, char* argv[]) {
  namespace po = boost::program_options;

  // define the general command line arguments
  command_line_arguments args;
  po::options_description arguments_description("Available options");
  std::size_t seed{};
  double time_limit_sec{};
  double observer_sec{};
  arguments_description.add_options()("help,h", "print this help message");
  arguments_description.add_options()("protein,p",
                                      po::value(&args.protein_path)->default_value(args.protein_path),
                                      "Path to the protein file (in PDB)");
  arguments_description.add_options()("ligand,l",
                                      po::value(&args.ligand_path)->default_value(args.ligand_path),
                                      "Path to the ligands file (in MOL2)");
  arguments_description.add_options()(
      "use",
      po::value<std::vector<std::string>>(&args.device_confs)->multitoken()->composing(),
      "Map each implementation to the device");
  arguments_description.add_options()(
      "time_limit_sec",
      po::value(&time_limit_sec),
      "Optional benchmark time limit in seconds; when reached, pending input ligands are discarded");
  arguments_description.add_options()(
      "observer",
      po::value(&observer_sec),
      "Optional throughput observer interval in seconds");
  // define the knobs command line arguments
  po::options_description knobs_description("Virtual Screening Knobs");
  knobs_description.add_options()(
      "population",
      po::value(&args.knobs.population_number)->default_value(args.knobs.population_number),
      "Number of individual(s) in the GA population");
  knobs_description.add_options()(
      "generations",
      po::value(&args.knobs.num_generations)->default_value(args.knobs.num_generations),
      "Number of generations that GA simulates");
  knobs_description.add_options()(
      "tournament_len",
      po::value(&args.knobs.tournament_length)->default_value(args.knobs.tournament_length),
      "Number of classes to select a parent in GA");
  knobs_description.add_options()(
      "mutation",
      po::value(&args.knobs.mutation_prob)->default_value(args.knobs.mutation_prob),
      "Probability of a mutation to happen during GA");
  knobs_description.add_options()("seed", po::value(&seed), "Seed for random values generators");
  knobs_description.add_options()(
      "tokens",
      po::value(&args.knobs.max_tbb_tokens)->default_value(args.knobs.max_tbb_tokens),
      "Max number of tokens in the TBB pipeline");
  knobs_description.add_options()(
      "bytes_per_token",
      po::value(&args.knobs.max_bytes_per_token)->default_value(args.knobs.max_bytes_per_token),
      "Max number of bytes per token in the TBB pipeline");
  knobs_description.add_options()(
      "queue_size",
      po::value(&args.knobs.max_tbb_queue_size)->default_value(args.knobs.max_tbb_queue_size),
      "Max number of ligands buffered in the TBB input/output queues");
  // parse them
  po::options_description all("Allowed Options");
  all.add(arguments_description).add(knobs_description);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);

  // handle the help message
  if (vm.count("help") > 0) {
    std::cout << "This application reads ligands from --ligand/-l and prints one score per output line."
              << std::endl;
    std::cout << std::endl;
    std::cout << "USAGE: " << argv[0] << " --protein|-p " << args.protein_path << " --ligand|-l "
              << args.ligand_path << " --use " << use_cpu_conf << " [MORE_CONFIGS...] [KNOBS] " << std::endl;
    std::cout << std::endl;
    std::cout << arguments_description << std::endl;
    std::cout << std::endl;
    std::cout << knobs_description << std::endl;
    std::cout << std::endl;
    std::cout << "The use flag accepts one or more configurations that describe which implementation" << std::endl
              << "should run on which hardware." << std::endl
              << "It has the following grammar: " << std::endl
              << "  --use CONFIGURATION [CONFIGURATION ...]" << std::endl
              << "  CONFIGURATION  -> IMPLEMENTATION:DEVICE:IDS[:WORKERS][:MEMORY_BYTES]" << std::endl
              << "  IMPLEMENTATION -> backend token such as CPP, CUDA, HIP, SYCL, GH, XSIMD" << std::endl
              << "  DEVICE         -> CPU|GPU" << std::endl
              << "  IDS            -> GROUP[,GROUP]*" << std::endl
              << "  GROUP          -> <device_id>|<device_id>-<device_id>" << std::endl
              << "  WORKERS        -> number of workers per GPU device (optional)" << std::endl
              << "  MEMORY_BYTES   -> per-device bucket memory budget in bytes (optional)" << std::endl
              << "The <device_id> number is directly related to the device id, while the option" << std::endl
              << "<device_id>-<device_id> can be used to specify a range" << std::endl;
    exit(EXIT_SUCCESS);
  }

  // make sure that the arguments make sense before returning them
  po::notify(vm);
  if (vm.count("seed")) {
    args.knobs.seed = std::optional<size_t>{seed};
  }
  if (vm.count("time_limit_sec")) {
    args.time_limit_sec = std::optional<double>{time_limit_sec};
  }
  if (vm.count("observer")) {
    args.observer = std::optional<double>{observer_sec};
  }
  return args;
}
