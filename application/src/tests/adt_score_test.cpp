#include "mudock/cpp_implementation/vectorization.hpp"

#include <boost/program_options.hpp>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <mudock/chem/autodock_ligand_types.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/grid/grid_map.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/fragments.hpp>
#include <mudock/mudock.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tests/autogrid.hpp>

template<class T>
inline T round3dp(const T x) {
  return ((std::floor((x) * 1000.0 + 0.5)) / 1000.0);
}

struct dpf_tokens {
  static constexpr auto MAP_TOKEN    = "map ";
  static constexpr auto ELEC_TOKEN   = "elecmap";
  static constexpr auto DESOLV_TOKEN = "dsolvmap";
  static constexpr auto MOVE_TOKEN   = "move";
};

static constexpr auto SCORE_TOKEN       = "AUTODOCK_SCORE";
static constexpr auto ERROR_SCORE_TOKEN = "AUTODOCK_ERROR_CORRECTION";

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

  mudock::info("Reading and parsing DPF file ", dpf_path, " ...");
  const auto desc = read_from_stream(std::ifstream(dpf_path));
  std::stringstream desc_s{desc};

  std::string line;
  std::unique_ptr<mudock::grid_map> electrostatic_map, desolvation_map;
  std::vector<std::unique_ptr<mudock::grid_atom_map>> grid_maps_unique;
  auto ligand = std::make_unique<mudock::static_molecule>();
  mudock::fp_type adt_score, adt_error_score;
  // TODO check if you can get rid of this and use mudock autogrid maps
  while (std::getline(desc_s, line)) {
    // Skip empty lines
    if (line.empty())
      continue;

    if (line.find(dpf_tokens::MAP_TOKEN) != std::string::npos) {
      std::stringstream ss{line};
      std::string map_path, _;
      ss >> _ >> map_path;
      if (line.find(dpf_tokens::ELEC_TOKEN) != std::string::npos) {
        electrostatic_map = std::make_unique<mudock::grid_map>(load_autogrid_map(map_path));
      } else if (line.find(dpf_tokens::DESOLV_TOKEN) != std::string::npos) {
        desolvation_map = std::make_unique<mudock::grid_map>(load_autogrid_map(map_path));
      } else {
        size_t last_dot_pos        = map_path.rfind('.');
        size_t second_last_dot_pos = map_path.rfind('.', last_dot_pos - 1);
        const auto symbol = map_path.substr(second_last_dot_pos + 1, last_dot_pos - second_last_dot_pos - 1);
        const auto map_type = mudock::autodock_type_from_ligand(mudock::parse_map_symbol(symbol));
        grid_maps_unique.emplace_back(
            std::make_unique<mudock::grid_atom_map>(map_type, load_autogrid_map(map_path)));
      }
    } else if (line.find(dpf_tokens::MOVE_TOKEN) != std::string::npos) {
      std::stringstream ss{line};
      std::string ligand_path, _;
      ss >> _ >> ligand_path;

      mudock::parse(*ligand, ligand_path);
      const auto ob_mol = mudock::parser(ligand_path);
      convert<mudock::pdbqt_rotate_check>(*ligand, ob_mol);
      mudock::apply_autodock_forcefield_pdbqt(*ligand, ligand_path);
    } else if (line.find(SCORE_TOKEN) != std::string::npos) {
      std::stringstream ss{line};
      std::string _;
      ss >> _ >> _ >> adt_score;
    } else if (line.find(ERROR_SCORE_TOKEN) != std::string::npos) {
      std::stringstream ss{line};
      std::string _;
      ss >> _ >> _ >> adt_error_score;
    }
  }

  std::sort(
      grid_maps_unique.begin(),
      grid_maps_unique.end(),
      [](const std::unique_ptr<mudock::grid_atom_map>& a, const std::unique_ptr<mudock::grid_atom_map>& b) {
        return static_cast<int>(autodock_ligand_from_type(a->get_atom_type())) <
               static_cast<int>(autodock_ligand_from_type(b->get_atom_type())); // descending order
      });

  std::vector<mudock::grid_atom_map> grid_maps;
  for (auto& grid: grid_maps_unique) grid_maps.push_back(*grid.release());
  const auto grid_atom_maps = std::make_shared<const mudock::grid_atom_mapper>(grid_maps);

  const auto num_atoms    = ligand->num_atoms();
  const int atom_map_size = grid_atom_maps.get()->get_single_map_size();

  std::vector<int> map_ligand_offsets;
  map_ligand_offsets.resize(num_atoms);
  for (int i = 0; i < num_atoms; i++)
    map_ligand_offsets[i] =
        static_cast<int>(autodock_ligand_from_type(ligand->autodock_type(i))) * atom_map_size;

  std::vector<int> non_bond_list_a1, non_bond_list_a2;
  std::vector<mudock::fp_type> cA_v, cB_v;
  std::vector<int> xB_v;
  mudock::non_bond_list(*ligand, non_bond_list_a1, non_bond_list_a2);
  mudock::precompute_lennard_jones(non_bond_list_a1.size(),
                                   cA_v,
                                   cB_v,
                                   xB_v,
                                   *ligand,
                                   non_bond_list_a1,
                                   non_bond_list_a2);

  mudock::info("Computing energy ...");
  const auto energy = mudock::calc_energy<mudock::cpu_vectorization::AUTO>(
      ligand->get_x().data(),
      ligand->get_y().data(),
      ligand->get_z().data(),
      ligand->get_vol().data(),
      ligand->get_solpar().data(),
      ligand->get_charge().data(),
      map_ligand_offsets.data(),
      num_atoms,
      ligand->num_rotamers(),
      non_bond_list_a1.size(),
      non_bond_list_a1.data(),
      non_bond_list_a2.data(),
      cA_v.data(),
      cB_v.data(),
      xB_v.data(),
      electrostatic_map.get()->minimum_coord.get_array().data(),
      electrostatic_map.get()->maximum_coord.get_array().data(),
      electrostatic_map.get()->center.get_array().data(),
      electrostatic_map.get()->index.size_x(),
      electrostatic_map.get()->index.size_xy(),
      grid_atom_maps->get_fused_maps().data(),
      electrostatic_map.get()->data(),
      desolvation_map.get()->data());

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
