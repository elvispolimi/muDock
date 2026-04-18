#include <boost/program_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <mudock/compute/safe_stack.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/format/writer.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/splitter.hpp>
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

  const auto in_file_format  = mudock::parse_supported_format(input_file);
  const auto out_file_format = mudock::parse_supported_format(output_file);
  constexpr_switch<0, mudock::get_num_supported_format(), 1>(
      [&](const auto format_index_in) {
        constexpr mudock::supported_format in_format =
            static_cast<mudock::supported_format>(format_index_in());
        auto input_text = read_from_stream(std::ifstream(input_file));
        mudock::splitter<mudock::type_of_format<in_format>> split;
        auto ligands_description = split(std::move(input_text));
        if (auto remainder = split.flush(); !remainder.empty()) {
          ligands_description.emplace_back(std::move(remainder));
        }
        // parse the input ligands and put them in a stack that we can compute
        mudock::info("Parsing ", ligands_description.size(), " compound(s) ...");
        std::size_t skipped_compounds = 0;
        constexpr_switch<0, mudock::get_num_supported_format(), 1>(
            [&](const auto format_index_out) {
              constexpr mudock::supported_format out_format =
                  static_cast<mudock::supported_format>(format_index_out());

              { std::ofstream ofs(output_file, std::ios::trunc); }
              std::ofstream ofs(output_file, std::ios::out | std::ios::app);

              for (const auto& description: ligands_description) {
                try {
                  mudock::writer<out_format, mudock::dynamic_molecule>(
                      mudock::parser<in_format, mudock::dynamic_molecule>(description),
                      ofs);
                } catch (...) {
                  ++skipped_compounds;
                }
              }
            },
            out_file_format);
        if (skipped_compounds > 0) {
          mudock::error("Skipped ", skipped_compounds, " compound(s) due to parse or write errors.");
        }
      },
      in_file_format);

  mudock::info("Converted into ", output_file);

  return EXIT_SUCCESS;
}
