#include "command_line_args.hpp"

#include <iostream>
#include <mudock/mudock.hpp>
#include <stdexcept>
#include <string>

int main(int argc, char* argv[]) {
  const auto args = parse_command_line_arguments(argc, argv);

  // find all the ligands description from the standard input
  auto input_text = std::string{std::istreambuf_iterator<std::string::value_type>{std::cin},
                                std::istreambuf_iterator<std::string::value_type>{}};
  mudock::splitter<mudock::mol2> split;
  auto ligands_description = split(std::move(input_text));
  ligands_description.emplace_back(split.flush());

  // parse the description to populate the actual data structures
  mudock::mol2 mol2;
  auto ligands = std::vector<mudock::static_molecule>{ligands_description.size()};
  for (std::size_t i{0}; i < ligands.size(); ++i) {
    try {
      mol2.parse(ligands[i], ligands_description[i]);
    } catch (const std::runtime_error& e) {
      std::cerr << "Unable to parse the ligand with index " << i << ", due to: " << e.what() << std::endl;
    }
  }

  // print the name of the molecules
  for (const auto& ligand: ligands) {
    const auto name = ligand.properties.get(mudock::property_type::NAME);
    std::cout << "Read ligand " << name << std::endl;
  }

  return EXIT_SUCCESS;
}
