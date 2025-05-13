#include <boost/program_options.hpp>
#include <command_line_args.hpp>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <mudock/chem/ligand_maps.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/grid/grid_map.hpp>
#include <mudock/knobs.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/fragments.hpp>
#include <mudock/mudock.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <sstream>
#include <string>
#include <tests/autogrid.hpp>

int main(int argc, char* argv[]) {
  namespace po                          = boost::program_options;
  std::filesystem::path protein_path    = std::filesystem::path{"protein.pdb"};
  std::filesystem::path ligand_path     = std::filesystem::path{"ligand.mol2"};
  const auto knobs                      = mudock::knobs{1, 1, 0, 0, 123};
  std::vector<std::string> device_confs = {std::string{use_cpu_conf}};

  po::options_description arguments_description("Available options");
  arguments_description.add_options()("help", "print this help message");
  arguments_description.add_options()("protein,p",
                                      po::value(&protein_path)->default_value(protein_path),
                                      "Path to the protein file");
  arguments_description.add_options()("ligand,l",
                                      po::value(&ligand_path)->default_value(ligand_path),
                                      "Path to the ligand file");
  arguments_description.add_options()(
      "use",
      po::value<std::vector<std::string>>(&device_confs)->multitoken()->composing(),
      "Map each implementation to the device");

  // parse them
  po::options_description all("Allowed Options");
  all.add(arguments_description);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);

  po::notify(vm);

  mudock::info("Reading and parsing protein ", protein_path, " ...");
  auto protein_ptr = std::make_shared<mudock::dynamic_molecule>();
  auto& protein    = *protein_ptr;
  parse(protein, protein_path);

  mudock::apply_autodock_forcefield(protein);
  auto grid_atom_maps    = std::make_shared<const mudock::grid_atom_mapper>(generate_atom_grid_maps(protein));
  auto electrostatic_map = std::make_shared<const mudock::grid_map>(generate_electrostatic_grid_map(protein));
  auto desolvation_map   = std::make_shared<const mudock::grid_map>(generate_desolvation_grid_map(protein));

  mudock::info("Reading and parsing ligand ", ligand_path, " ...");
  auto ligand = std::make_unique<mudock::static_molecule>();
  mudock::parse(*ligand, ligand_path);

  mudock::apply_autodock_forcefield(*ligand);

  mudock::info("Generating score reference ...");
  auto output_queue = std::make_shared<mudock::safe_stack<mudock::static_molecule>>();
  auto input_queue  = std::make_shared<mudock::safe_stack<mudock::static_molecule>>();
  input_queue->enqueue(std::make_unique<mudock::static_molecule>(mudock::static_molecule(*ligand)));
  {
    auto threadpool = mudock::threadpool();
    mudock::manage_cpp({std::string{use_cpu_conf}},
                       threadpool,
                       grid_atom_maps,
                       electrostatic_map,
                       desolvation_map,
                       knobs,
                       input_queue,
                       output_queue);
  }
  auto ligand_out = output_queue->dequeue();
  std::stringstream ss(ligand_out->properties.get(mudock::property_type::SCORE));
  mudock::fp_type reference_score;
  ss >> reference_score;

  for (auto& conf: device_confs) {
    mudock::info("Comparing reference with ", conf, " ...");
    input_queue->enqueue(std::make_unique<mudock::static_molecule>(mudock::static_molecule(*ligand)));
    {
      auto threadpool = mudock::threadpool();
      mudock::manage_cpp({conf},
                         threadpool,
                         grid_atom_maps,
                         electrostatic_map,
                         desolvation_map,
                         knobs,
                         input_queue,
                         output_queue);
      mudock::manage_cuda({conf},
                          threadpool,
                          knobs,
                          grid_atom_maps,
                          electrostatic_map,
                          desolvation_map,
                          input_queue,
                          output_queue);
      mudock::manage_hip({conf},
                         threadpool,
                         knobs,
                         grid_atom_maps,
                         electrostatic_map,
                         desolvation_map,
                         input_queue,
                         output_queue);
      mudock::manage_sycl({conf},
                          threadpool,
                          knobs,
                          grid_atom_maps,
                          electrostatic_map,
                          desolvation_map,
                          input_queue,
                          output_queue);
      mudock::manage_omp({conf},
                         threadpool,
                         knobs,
                         grid_atom_maps,
                         electrostatic_map,
                         desolvation_map,
                         input_queue,
                         output_queue);
    }
    ligand_out = output_queue->dequeue();
    std::stringstream sss(ligand_out->properties.get(mudock::property_type::SCORE));
    mudock::fp_type score;
    sss >> score;
    const auto diff  = std::abs(score - reference_score);
    const auto error = std::max(reference_score * mudock::fp_type{0.001}, mudock::fp_type{0.1});
    if (diff > error) {
      mudock::error(std::format("Difference betweem scores of {} on {} ( CPU {} vs {} {})",
                                ligand_path.string(),
                                protein_path.string(),
                                reference_score,
                                conf,
                                score));
      throw std::runtime_error("Error in score");
    }
  }

  mudock::info(
      std::format("Succesfully verified the score of {} on {}", ligand_path.string(), protein_path.string()));
  return EXIT_SUCCESS;
}
