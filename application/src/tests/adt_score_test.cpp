#include "autogrid.hpp"

#include <boost/program_options.hpp>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/adt_score.hpp>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/log.hpp>
#include <mudock/mudock.hpp>
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
  mudock::autodock_ligand adt_ligand{ligand};
  const auto adt_score       = load_autodock_score(dpf_path);
  const auto adt_error_score = load_autodock_error_score(dpf_path);

  const auto num_atoms    = ligand.num_atoms();
  const auto num_rotamers = ligand.num_rotamers();
  adt_ligand.update_offsets(static_cast<int>(adt_grid.get_map_flat_size()));

  mudock::info("Computing energy ...");
  std::vector<int> num_atoms_b{num_atoms};
  std::vector<int> num_rotamers_b{static_cast<int>(num_rotamers)};
  std::vector<int> num_nonbonds_b{0, static_cast<int>(adt_ligand.non_bond_size())};
  std::vector<mudock::fp_type> scores_b{0};
  auto q = std::make_shared<mudock::queue_cpp>(0, mudock::device_type::CPU);
  mudock::adt_score_kernel<mudock::queue_cpp> adt_kernel{1,
                                                         1,
                                                         num_atoms,
                                                         num_atoms_b.data(),
                                                         num_rotamers_b.data(),
                                                         num_nonbonds_b.data(),
                                                         ligand.x(),
                                                         ligand.y(),
                                                         ligand.z(),
                                                         adt_ligand.vol(),
                                                         adt_ligand.solpar(),
                                                         ligand.charge(),
                                                         adt_ligand.atom_map_offsets(),
                                                         adt_ligand.non_bond_A(),
                                                         adt_ligand.non_bond_B(),
                                                         adt_ligand.non_bond_cA(),
                                                         adt_ligand.non_bond_cB(),
                                                         adt_ligand.non_bond_xB(),
                                                         adt_grid.get_maps_pointer(),
                                                         adt_grid.get_min_p(),
                                                         adt_grid.get_max_p(),
                                                         adt_grid.get_center_p(),
                                                         static_cast<int>(adt_grid.get_size_x()),
                                                         static_cast<int>(adt_grid.get_size_xy()),
                                                         static_cast<int>(adt_grid.get_size_xyz()),
                                                         scores_b.data(),
                                                         q};

  adt_kernel();
  const mudock::fp_type energy = scores_b[0];
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
