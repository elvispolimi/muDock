#include "command_line_args.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <mudock/mudock.hpp>
#include <stdexcept>
#include <string>

// utility function that reads the whole content of a stream
template<class stream_type>
inline auto read_from_stream(stream_type&& in) {
  assert(in.good());
  return std::string{std::istreambuf_iterator<std::string::value_type>{in},
                     std::istreambuf_iterator<std::string::value_type>{}};
}

int main(int argc, char* argv[]) {
  const auto args = parse_command_line_arguments(argc, argv);

  // read and parse the target protein
  auto protein_ptr               = std::make_shared<mudock::dynamic_molecule>();
  auto& protein                  = *protein_ptr;
  auto pdb                       = mudock::pdb{};
  const auto protein_description = read_from_stream(std::ifstream(args.protein_path));
  pdb.parse(protein, protein_description);
  mudock::apply_autodock_forcefield(protein);

  // read  all the ligands description from the standard input and split them
  auto input_text = read_from_stream(std::cin);
  mudock::splitter<mudock::mol2> split;
  auto ligands_description = split(std::move(input_text));
  ligands_description.emplace_back(split.flush());

  // parse the input ligands and put them in a stack that we can compute
  mudock::mol2 mol2;
  auto input_queue = mudock::safe_stack<mudock::static_molecule>{};
  for (const auto& description: ligands_description) {
    try {
      auto ligand = std::make_unique<mudock::static_molecule>();
      mol2.parse(*ligand, description);
      mudock::apply_autodock_forcefield(*ligand);
      input_queue.enqueue(std::move(ligand));
    } catch (const std::runtime_error& e) {
      std::cerr << "Unable to parse the following ligand: " << std::endl;
      std::cerr << description << std::endl;
      std::cerr << "Due to: " << e.what() << std::endl;
    }
  }

  // compute all the ligands using a single cpp implementation
  auto output_queue = mudock::safe_stack<mudock::static_molecule>{};
  auto worker       = mudock::cpp_worker(protein_ptr, input_queue, output_queue, 0);
  worker.main();

  return EXIT_SUCCESS;
}
