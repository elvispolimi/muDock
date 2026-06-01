#include <boost/program_options.hpp>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <format>
#include <iostream>
#include <command_line_args.hpp>
#include <mudock/chem/autodock_layer.hpp>
#include <mudock/compute/manager.hpp>
#include <mudock/compute/pipeline.hpp>
#include <mudock/cpp_implementation/vinardo_score_kernel_cpp.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/containers.hpp>
#include <stdexcept>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace {

std::unordered_map<std::string, mudock::fp_type>
compute_vinardo_affinities(const std::filesystem::path& receptor_path,
                           const std::vector<std::filesystem::path>& ligand_paths) {
  mudock::info("Reading and parsing protein ", receptor_path, " ...");
  auto protein =
      std::make_shared<mudock::dynamic_molecule>(mudock::parser<mudock::dynamic_molecule>(receptor_path));

  mudock::knobs conf;
  conf.population_number = 1;
  // Instantiate the scoring pipeline and the queues
  // as done in the adt test
  mudock::vinardo_score_pipeline pipeline{protein};
  auto output_queue = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
  auto input_queue  = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
  input_queue->initialize(ligand_paths.size());
  output_queue->initialize(ligand_paths.size());

  [[maybe_unused]] mudock::autodock_dynamic_layer protein_autodock{
      *protein,
      [receptor_path](mudock::autodock_dynamic_layer& layer) {
        mudock::apply_autodock_forcefield_pdbqt(layer, receptor_path);
      }};

  for (const auto& ligand_path: ligand_paths) {
    mudock::info("Reading and parsing ligand ", ligand_path, " ...");
    auto ligand = std::make_unique<mudock::static_molecule>(
        mudock::parser<mudock::static_molecule>(ligand_path, &mudock::pdbqt_rotate_check));

    [[maybe_unused]] mudock::autodock_static_layer ligand_autodock{
        *ligand,
        [ligand_path](mudock::autodock_static_layer& layer) {
          mudock::apply_autodock_forcefield_pdbqt(layer, ligand_path);
        }};

    if (!ligand->pdbqt_ligand_data.valid) {
      throw std::runtime_error("Missing PDBQT ligand data");
    }

    ligand->properties.assign(mudock::property_type::NAME, ligand_path.string());
    input_queue->enqueue(ligand);
  }
  input_queue->send_terminate_signal();
  {
    auto threadpool = mudock::threadpool();
    mudock::manager({std::string{use_cpu_conf}}, threadpool, conf, input_queue, output_queue, pipeline);
  }
  output_queue->send_terminate_signal();

  std::unordered_map<std::string, mudock::fp_type> scores;
  for (auto ligand_out = output_queue->dequeue(); ligand_out; ligand_out = output_queue->dequeue()) {
    std::stringstream ss{ligand_out->properties.get(mudock::property_type::SCORE)};
    mudock::fp_type score{};
    ss >> score;
    scores.emplace(ligand_out->properties.get(mudock::property_type::NAME), score);
  }

  if (scores.size() != ligand_paths.size()) {
    throw std::runtime_error("Missing scored ligand");
  }

  return scores;
}

} // namespace

int main(int argc, char* argv[]) {
  namespace po = boost::program_options;

  std::filesystem::path receptor_path;
  std::vector<std::filesystem::path> ligand_paths;
  std::vector<std::filesystem::path> reference_paths;
  mudock::fp_type tolerance = mudock::fp_type{1e-3};

  po::options_description arguments_description("Available options");
  arguments_description.add_options()("help", "print this help message");
  arguments_description.add_options()("receptor,r",
                                      po::value(&receptor_path),
                                      "Path to the receptor PDBQT file");
  arguments_description.add_options()("ligand,l",
                                      po::value(&ligand_paths)->multitoken()->composing(),
                                      "Path to the ligand PDBQT file");
  arguments_description.add_options()(
      "reference", po::value(&reference_paths)->multitoken()->composing(), "Path to reference affinity");
  arguments_description.add_options()(
      "tolerance", po::value(&tolerance)->default_value(tolerance), "Absolute tolerance");

  po::options_description all("Allowed Options");
  all.add(arguments_description);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
  po::notify(vm);

  if (vm.contains("help")) {
    std::cout << all << '\n';
    return EXIT_SUCCESS;
  }

  if (receptor_path.empty() || ligand_paths.empty()) {
    throw std::runtime_error("Both receptor and ligand paths are required");
  }

  const auto mudock_affinities = compute_vinardo_affinities(receptor_path, ligand_paths);
  if (reference_paths.empty()) {
    for (const auto& ligand_path: ligand_paths) {
      mudock::info(std::format("Vinardo affinity for {}: {}",
                               ligand_path.string(),
                               mudock_affinities.at(ligand_path.string())));
    }
    return EXIT_SUCCESS;
  }

  if (reference_paths.size() != ligand_paths.size()) {
    throw std::runtime_error("Number of references must match number of ligands");
  }

  for (std::size_t i = 0; i < ligand_paths.size(); ++i) {
    const auto& ligand_path    = ligand_paths[i];
    const auto& reference_path = reference_paths[i];
    std::ifstream input{reference_path};
    if (!input) {
      throw std::runtime_error(std::format("Cannot open reference file {}", reference_path.string()));
    }

    mudock::fp_type smina_affinity{};
    input >> smina_affinity;
    const auto mudock_affinity = mudock_affinities.at(ligand_path.string());
    const auto diff            = mudock_affinity - smina_affinity;
    mudock::info(std::format("muDock={} smina={} diff={}", mudock_affinity, smina_affinity, diff));

    if (std::abs(diff) > tolerance) {
      mudock::error(std::format("{} exceeds tolerance {} with abs diff {}",
                                ligand_path.string(),
                                tolerance,
                                std::abs(diff)));
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
