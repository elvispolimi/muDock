#include "command_line_args.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <mudock/mudock.hpp>
#include <stdexcept>
#include <string>

int main(int argc, char* argv[]) {
  const auto args = parse_command_line_arguments(argc, argv);

  // read and parse the target protein
  auto protein      = mudock::dynamic_molecule{};
  auto pdb          = mudock::pdb{};
  auto protein_file = std::ifstream(args.protein_path);
  assert(protein_file.good());
  const auto protein_description =
      std::string{std::istreambuf_iterator<std::string::value_type>{protein_file},
                  std::istreambuf_iterator<std::string::value_type>{}};
  pdb.parse(protein, protein_description);

  // find all the ligands description from the standard input
  auto input_text = std::string{std::istreambuf_iterator<std::string::value_type>{std::cin},
                                std::istreambuf_iterator<std::string::value_type>{}};
  mudock::splitter<mudock::mol2> split;
  auto ligands_description = split(std::move(input_text));
  ligands_description.emplace_back(split.flush());

  // parse the description to populate the actual data structures
  mudock::mol2 mol2;
  auto i_queue = std::make_shared<mudock::squeue<mudock::static_molecule>>();
  for (std::size_t i{0}; i < ligands_description.size(); ++i) {
    try {
      std::unique_ptr<mudock::static_molecule> ligand = std::make_unique<mudock::static_molecule>();
      mol2.parse(*(ligand.get()), ligands_description[i]);
      mudock::apply_autodock_forcefield(*ligand.get());

      i_queue->enqueue(std::move(ligand));
    } catch (const std::runtime_error& e) {
      std::cerr << "Unable to parse the ligand with index " << i << ", due to: " << e.what() << std::endl;
    }
  }

  const mudock::device_conf d_c = mudock::parse_conf(args.device_conf);
  auto o_queue                  = std::make_shared<mudock::squeue<mudock::static_molecule>>();

  constexpr_for<static_cast<int>(mudock::kernel_type::CPP), static_cast<int>(mudock::kernel_type::NONE), 1>(
      [&](auto index) {
        constexpr mudock::kernel_type sel_l = static_cast<mudock::kernel_type>(index.value);
        if (sel_l == d_c.l_t) {
          // NOTE: The destructor will handle the termination
          mudock::thread_pool<sel_l> t_pool{i_queue, o_queue, d_c};
        }
      });

  std::unique_ptr<mudock::static_molecule> mol = o_queue->dequeue();
  while (mol.get() != nullptr) {
    const auto name  = mol->properties.get(mudock::property_type::NAME);
    const auto score = mol->properties.get(mudock::property_type::SCORE);
    std::cout << "Read ligand " << name << " with score " << score << std::endl;
    mol = o_queue->dequeue();
  }

  return EXIT_SUCCESS;
}
