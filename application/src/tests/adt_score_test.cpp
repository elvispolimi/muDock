#include "autogrid.hpp"

#include <boost/program_options.hpp>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <format>
#include <memory>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/adt_score.hpp>
#include <mudock/compute/buffer.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/log.hpp>
#include <mudock/mudock.hpp>
#include <mudock/knobs.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <stdexcept>
#include <string>

template<class T>
inline T round3dp(const T x) {
  return ((std::floor((x) * 1000.0 + 0.5)) / 1000.0);
}

int main(int argc, char *argv[]) {
  namespace po                   = boost::program_options;
  std::filesystem::path dpf_path = std::filesystem::path{"rec.dpf"};

  po::options_description arguments_description("Available options");
  arguments_description.add_options()("help", "print this help message");
  arguments_description.add_options()("dpf",
                                      po::value(&dpf_path)->default_value(dpf_path),
                                      "Path to the autodock DPF file");

  // parse them
  po::options_description all("Allowed Options");
  all.add(arguments_description);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);

  po::notify(vm);

  mudock::autodock_grid adt_grid = load_autogrid_map_dpf(dpf_path);

  const auto ligand_path = get_ligand_path(dpf_path);
  const auto ligand_format = mudock::parse_supported_format(std::filesystem::path{ligand_path});
  const auto rotate_check  = ligand_format == mudock::supported_format::PDBQT ? &mudock::pdbqt_rotate_check
                                                                               : &mudock::ob_rotate_check;
  mudock::static_molecule ligand = mudock::parser<mudock::static_molecule>(ligand_path, rotate_check);
  const auto adt_score       = load_autodock_score(dpf_path);
  const auto adt_error_score = load_autodock_error_score(dpf_path);

  mudock::info("Computing energy ...");
  mudock::knobs configuration{};
  configuration.population_number = 1;
  auto scratch = std::make_shared<mudock::scratchpad<mudock::queue_cpp>>(
      configuration, 0, mudock::device_type::CPU);
  auto device_scratch = std::make_shared<mudock::scratchpad<mudock::queue_cpp>>(
      configuration, 0, mudock::device_type::CPU);

  auto& prot_min = device_scratch->get<mudock::buffer_data_type::PROT_MIN>();
  auto& prot_max = device_scratch->get<mudock::buffer_data_type::PROT_MAX>();
  auto& prot_center = device_scratch->get<mudock::buffer_data_type::PROT_CENTER>();
  auto& size_x = device_scratch->get<mudock::buffer_data_type::PROT_SIZE_X>();
  auto& size_xy = device_scratch->get<mudock::buffer_data_type::PROT_SIZE_XY>();
  auto& size_xyz = device_scratch->get<mudock::buffer_data_type::PROT_SIZE_XYZ>();
  auto& maps = device_scratch->get<mudock::buffer_data_type::PROT_GRID_MAPS>();
  prot_min.alloc(3);
  prot_max.alloc(3);
  prot_center.alloc(3);
  size_x.alloc(1);
  size_xy.alloc(1);
  size_xyz.alloc(1);
  maps.alloc(adt_grid.get_map_flat_size() * mudock::num_autodock_grids());
  std::memcpy(prot_min(), adt_grid.get_min_p(), 3 * sizeof(mudock::fp_type));
  std::memcpy(prot_max(), adt_grid.get_max_p(), 3 * sizeof(mudock::fp_type));
  std::memcpy(prot_center(), adt_grid.get_center_p(), 3 * sizeof(mudock::fp_type));
  size_x()[0] = static_cast<int>(adt_grid.get_size_x());
  size_xy()[0] = static_cast<int>(adt_grid.get_size_xy());
  size_xyz()[0] = static_cast<int>(adt_grid.get_size_xyz());
  std::memcpy(maps(), adt_grid.get_maps_pointer(),
              adt_grid.get_map_flat_size() * mudock::num_autodock_grids() * sizeof(mudock::fp_type));

  mudock::dynamic_molecule protein;
  protein.resize(1, 0);
  protein.x(0) = adt_grid.get_center_p()[0];
  protein.y(0) = adt_grid.get_center_p()[1];
  protein.z(0) = adt_grid.get_center_p()[2];
  mudock::adt_score<mudock::queue_cpp> scorer{
      scratch,
      device_scratch,
      protein,
      mudock::ligand_placement{mudock::ligand_placement_mode::preserve, std::nullopt}};
  mudock::batch<mudock::static_molecule> batch;
  batch.num_ligands = 1;
  batch.num_atoms = ligand.num_atoms();
  batch.batch_max_atoms = ligand.num_atoms();
  batch.molecules[0] = std::make_unique<mudock::static_molecule>(ligand);
  scorer.prepare(batch);
  scorer();
  const mudock::fp_type energy = scratch->get<mudock::buffer_data_type::SCORES>()()[0];
  // High tolerance due to the precomputation done in autodock, refers to intnbtable.cc
  if (std::abs(energy - adt_score + adt_error_score) > static_cast<mudock::fp_type>(0.1)) {
    // High tolerance due to the precomputation done in autodock, refers to intnbtable.cc
    if (std::abs(energy - adt_score + adt_error_score) > static_cast<mudock::fp_type>(0.1)) {
      mudock::error(std::format("Difference betweem scores of {} ({} vs {})",
                                dpf_path.string(),
                                adt_score - adt_error_score,
                                energy));
      throw std::runtime_error("Error in score");
    }
  }

  mudock::info(std::format("Succesfully verified the score in {}", dpf_path.string()));
  return EXIT_SUCCESS;
}
