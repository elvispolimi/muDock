#include "command_line_args.hpp"
#include "mudock/molecule.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/likwid_utils.hpp>
#include <mudock/mudock.hpp>
#include <string>

int main(int argc, char* argv[]) {
  const auto args = parse_command_line_arguments(argc, argv);

  // read and parse the target protein
  mudock::info("Reading and parsing protein ", args.protein_path, " ...");
  auto protein = std::make_shared<mudock::autodock_protein>();
  parse(*protein, args.protein_path);
  // TODO prepare protein

  // mudock::apply_autodock_forcefield(protein);
  // mudock::autodock_protein protein_adt = mudock::make_autodock_protein(protein);

  // read  all the ligands description from the standard input and split them
  mudock::info("Reading ligand ", args.ligand_path, " ...");
  const auto in_format = mudock::parse_supported_format(args.ligand_path);
  auto input_queue     = std::make_shared<mudock::safe_stack<mudock::autodock_ligand>>();
  constexpr_switch<0, mudock::get_num_supported_format(), 1>(
      [&](const auto format_index) {
        const auto format = static_cast<mudock::supported_format>(format_index());
        auto input_text   = read_from_stream(std::ifstream(args.ligand_path));
        mudock::splitter<mudock::type_of_format<static_cast<mudock::supported_format>(format_index())>> split;
        auto ligands_description = split(std::move(input_text));
        ligands_description.emplace_back(split.flush());

        // parse the input ligands and put them in a stack that we can compute
        mudock::info("Parsing ", ligands_description.size(), " ligand(s) ...");
        if constexpr (format == mudock::supported_format::MOL2X) {
#pragma omp parallel for shared(input_queue)
          for (const auto& description: ligands_description) {
            auto ligand = std::make_unique<mudock::autodock_ligand>();
            mudock::parse<format>(*ligand, description);
            ligand->prepare();
            // mudock::apply_autodock_forcefield(*ligand);
            input_queue->enqueue(ligand);
          }
        } else {
          for (const auto& description: ligands_description) {
            try {
              mudock::static_molecule ligand;
              mudock::parse<format>(ligand, description);
              // mudock::apply_autodock_forcefield(ligand);
              input_queue->enqueue(std::make_unique<mudock::autodock_ligand>(ligand));
            } catch (...) {}
          }
        }
      },
      in_format);

  // compute all the ligands according to the input configuration
  mudock::info("Virtual screening the ligands ...");
  LIKWID_MARKER_INIT;

  auto output_queue = std::make_shared<mudock::safe_stack<mudock::static_molecule>>();
  {
    auto threadpool = mudock::threadpool();
    mudock::manage_cpp(args.device_confs, threadpool, *protein, args.knobs, input_queue, output_queue);
    // mudock::manage_cuda(args.device_confs, threadpool, args.knobs, protein_adt, input_queue, output_queue);
    // mudock::manage_hip(args.device_confs, threadpool, args.knobs, protein_adt, input_queue, output_queue);
    // mudock::manage_sycl(args.device_confs, threadpool, args.knobs, protein_adt, input_queue, output_queue);
    // mudock::manage_omp(args.device_confs, threadpool, args.knobs, protein_adt, input_queue, output_queue);
    mudock::info("All workers have been created!");
  } // when we exit from this block the computation is complete

  LIKWID_MARKER_CLOSE;

  // after the computation it will be nice to print the score of all the molecules
  mudock::info("Printing the scores ...");
  for (auto ligand = output_queue->dequeue(); ligand; ligand = output_queue->dequeue()) {
    std::cout << ligand->properties.get(mudock::property_type::NAME) << " "
              << ligand->properties.get(mudock::property_type::SCORE) << std::endl;
  }

  // if we reach this statement we completed successfully the run
  mudock::info("All Done!");
  return EXIT_SUCCESS;
}
