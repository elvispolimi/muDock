#include <boost/program_options.hpp>
#include <cstdlib>
#include <filesystem>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/fragments.hpp>
#include <mudock/mudock.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <stdexcept>
#include <string>
#include <tests/autogrid.hpp>

template<class T>
inline T round3dp(const T x) {
  return ((std::floor((x) * 1000.0 + 0.5)) / 1000.0);
}

int main(int argc, char* argv[]) {
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

  mudock::autodock_protein protein = load_autogrid_map_dpf(dpf_path);
  auto ligand                      = load_autogrid_ligand(dpf_path);
  const auto adt_score             = load_autodock_score(dpf_path);
  const auto adt_error_score       = load_autodock_error_score(dpf_path);

  const auto num_atoms    = ligand.num_atoms();
  const auto num_rotamers = ligand.num_rotamers();
  ligand.update_offsets(protein);

  mudock::info("Computing energy ...");
  const auto energy = mudock::calc_energy<mudock::cpu_vectorization::AUTO>(ligand.x(),
                                                                           ligand.y(),
                                                                           ligand.z(),
                                                                           ligand.vol(),
                                                                           ligand.solpar(),
                                                                           ligand.charge(),
                                                                           ligand.atom_map_offsets(),
                                                                           num_atoms,
                                                                           num_rotamers,
                                                                           ligand.non_bond_size(),
                                                                           ligand.non_bond_A(),
                                                                           ligand.non_bond_B(),
                                                                           ligand.non_bond_cA(),
                                                                           ligand.non_bond_cB(),
                                                                           ligand.non_bond_xB(),
                                                                           protein.get_min_p(),
                                                                           protein.get_max_p(),
                                                                           protein.get_center_p(),
                                                                           protein.get_size_x(),
                                                                           protein.get_size_xy(),
                                                                           protein.get_size_xyz(),
                                                                           protein.get_maps_pointer());

  // High tolerance due to the precomputation done in autodock, refers to intnbtable.cc
  if (std::abs(energy - adt_score + adt_error_score) > mudock::fp_type{0.1}) {
    mudock::error(std::format("Difference betweem scores of {} ({} vs {})",
                              dpf_path.string(),
                              adt_score - adt_error_score,
                              energy));
    throw std::runtime_error("Error in score");
  }

  mudock::info(std::format("Succesfully verified the score in {}", dpf_path.string()));
  return EXIT_SUCCESS;
}
