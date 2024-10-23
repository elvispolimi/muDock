#include "command_line_args.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
// #include <mudock/cuda_implementation/map_textures.cuh>
#include <mudock/mudock.hpp>
#include <stdexcept>
#include <string>
#include <vector>

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
  mudock::info("Reading and parsing protein ", args.protein_path, " ...");
  auto protein_ptr               = std::make_shared<mudock::dynamic_molecule>();
  auto& protein                  = *protein_ptr;
  auto pdb                       = mudock::pdb{};
  const auto protein_description = read_from_stream(std::ifstream(args.protein_path));
  pdb.parse(protein, protein_description);
  ///////////////
  // std::ifstream infile("/home/gianmarco/ht_dataset/1a30_pocket.pdbqt");

  // std::string line;
  // size_t index{0};
  // while (std::getline(infile, line)) {
  //   // Only process lines that start with "ATOM" or "HETATM"
  //   if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
  //     // Extract atom type and charge
  //     protein.charge(index) = std::stod(line.substr(70, 6)); // Charge (column 71-76)

  //     index++;
  //   }
  // }
  ////////////////////

  mudock::apply_autodock_forcefield(protein);
  auto grid_atom_maps    = std::make_shared<const mudock::grid_atom_mapper>(generate_atom_grid_maps(protein));
  auto electrostatic_map = std::make_shared<const mudock::grid_map>(generate_electrostatic_grid_map(protein));
  auto desolvation_map   = std::make_shared<const mudock::grid_map>(generate_desolvation_grid_map(protein));
  // auto tp  = mudock::point3D{};
  // auto itp = mudock::index3D{1, 1, 1};
  // auto sv  = std::vector<mudock::grid_atom_map>{};
  // for (int i = 0; i < mudock::num_device_map_textures(); ++i)
  //   sv.emplace_back(mudock::autodock_type_from_map(static_cast<mudock::device_map_textures>(i)), protein);
  // auto grid_atom_maps    = std::make_shared<const mudock::grid_atom_mapper>(sv);
  // auto electrostatic_map = std::make_shared<const mudock::grid_map>(protein);
  // auto desolvation_map   = std::make_shared<const mudock::grid_map>(protein);

  // read  all the ligands description from the standard input and split them
  mudock::info("Reading ligands from the stdin ...");
  auto input_text = read_from_stream(std::cin);
  mudock::splitter<mudock::mol2> split;
  auto ligands_description = split(std::move(input_text));
  ligands_description.emplace_back(split.flush());

  // parse the input ligands and put them in a stack that we can compute
  mudock::info("Parsing ", ligands_description.size(), " ligand(s) ...");
  mudock::mol2 mol2;
  auto input_queue = std::make_shared<mudock::safe_stack<mudock::static_molecule>>();
  for (const auto& description: ligands_description) {
    try {
      auto ligand = std::make_unique<mudock::static_molecule>();
      mol2.parse(*ligand, description);
      ///////////////
      // infile = std::ifstream{"/home/gianmarco/ht_dataset/1a30_ligand.pdbqt"};

      // while (std::getline(infile, line)) {
      //   // Only process lines that start with "ATOM" or "HETATM"
      //   if (line.substr(0, 4) == "ATOM") {
      //     // Extract atom type and charge
      //     mudock::fp_type x = std::stod(line.substr(32, 6));
      //     mudock::fp_type y = std::stod(line.substr(40, 6));
      //     mudock::fp_type z = std::stod(line.substr(49, 6));
      //     for (index = 0; std::abs(ligand.get()->x(index) - x) > mudock::fp_type{0.001} ||
      //                     std::fabs(ligand.get()->y(index) - y) > mudock::fp_type{0.001} ||
      //                     std::fabs(ligand.get()->z(index) - z) > mudock::fp_type{0.001};
      //          ++index);
      //     ligand.get()->charge(index) = std::stod(line.substr(70, 6)); // Charge (column 71-76)

      //     index++;
      //   }
      // }
      ////////////////////
      mudock::apply_autodock_forcefield(*ligand);
      input_queue->enqueue(std::move(ligand));
    } catch (const std::exception& e) {
      std::cerr << "Unable to parse the following ligand: " << std::endl;
      std::cerr << description << std::endl;
      std::cerr << "Due to: " << e.what() << std::endl;
    }
  }

  // compute all the ligands according to the input configuration
  mudock::info("Virtual screening the ligands ...");
  auto output_queue = std::make_shared<mudock::safe_stack<mudock::static_molecule>>();
  {
    auto threadpool = mudock::threadpool();
    mudock::manage_cpp(args.device_conf,
                       threadpool,
                       grid_atom_maps,
                       electrostatic_map,
                       desolvation_map,
                       args.knobs,
                       input_queue,
                       output_queue);
    mudock::manage_cuda(args.device_conf,
                        threadpool,
                        args.knobs,
                        grid_atom_maps,
                        electrostatic_map,
                        desolvation_map,
                        input_queue,
                        output_queue);
    mudock::manage_hip(args.device_conf,
                       threadpool,
                       args.knobs,
                       grid_atom_maps,
                       electrostatic_map,
                       desolvation_map,
                       input_queue,
                       output_queue);
    mudock::manage_sycl(args.device_conf,
                        threadpool,
                        args.knobs,
                        grid_atom_maps,
                        electrostatic_map,
                        desolvation_map,
                        input_queue,
                        output_queue);
    mudock::manage_omp(args.device_conf,
                       threadpool,
                       args.knobs,
                       grid_atom_maps,
                       electrostatic_map,
                       desolvation_map,
                       input_queue,
                       output_queue);
    mudock::info("All workers have been created!");
  } // when we exit from this block the computation is complete

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
