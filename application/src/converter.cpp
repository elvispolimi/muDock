#include <boost/program_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/log.hpp>
#include <stdexcept>
#include <string>

int main(int argc, char* argv[]) {
  namespace po                      = boost::program_options;
  std::filesystem::path input_file  = std::filesystem::path{""};
  std::filesystem::path output_file = std::filesystem::path{""};

  po::options_description arguments_description("Argument descriptions");
  arguments_description.add_options()("help,h", "print this help message");
  arguments_description.add_options()("input,i", po::value(&input_file), "Path to the input file");
  arguments_description.add_options()("output,o", po::value(&output_file), "Path to the output file");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(arguments_description).run(), vm);

  po::notify(vm);
  if (vm.count("help") > 0) {
    std::cout << arguments_description << std::endl;
    return EXIT_SUCCESS;
  }
  if (!vm.count("input")) {
    throw std::runtime_error("Error: --input is required!");
  }
  if (!vm.count("output")) {
    throw std::runtime_error("Error: --output is required!");
  }

  mudock::info("Reading and parsing ", input_file, " ...");
  mudock::writer(mudock::parser(input_file), output_file);
  mudock::info("Converted into ", output_file);

  return EXIT_SUCCESS;
}
