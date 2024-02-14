#include "command_line_args.hpp"

#include <boost/program_options.hpp>
#include <iostream>

command_line_arguments parse_command_line_arguments(const int argc, char* argv[]) {
  namespace po = boost::program_options;

  // define the command lines
  command_line_arguments args;
  po::options_description arguments_description("Available options");
  arguments_description.add_options()("help", "print this help message");
  arguments_description.add_options()("protein",
                                      po::value(&args.protein_path)->default_value(args.protein_path),
                                      "Path to the protein file (in PDB)");

  // parse them
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(arguments_description).run(), vm);

  // handle the help message
  if (vm.count("help") > 0) {
    std::cout << "This application reads from the standard input a ligand library in mol2 format. It will"
              << std::endl;
    std::cout << "print on the standard output the best docked pose with the score as comment" << std::endl;
    std::cout << std::endl;
    std::cout << "USAGE: " << argv[0] << " --protein " << args.protein_path << " < \"/path/to/ligands.mol2\""
              << std::endl;
    std::cout << std::endl;
    std::cout << arguments_description << std::endl;
    std::cout << std::endl;
    exit(EXIT_SUCCESS);
  }

  // make sure that the arguments make sense before returning them
  po::notify(vm);
  return args;
}
