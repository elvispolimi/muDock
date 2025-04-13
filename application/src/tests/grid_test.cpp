#include <boost/program_options.hpp>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <mudock/chem/ligand_maps.hpp>
#include <mudock/format.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/grid/grid_map.hpp>
#include <mudock/log.hpp>
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

struct fld_tokens {
  static constexpr auto FILE_TOKEN     = "file=";
  static constexpr auto VARIABLE_TOKEN = "variable";
  static constexpr auto LABEL_TOKEN    = "label";
  static constexpr auto ELETRO_TOKEN   = "Electrostatics";
  static constexpr auto DESOLV_TOKEN   = "Desolvation";
};

int main(int argc, char* argv[]) {
  namespace po                     = boost::program_options;
  std::filesystem::path pdbqt_path = std::filesystem::path{"protein.pdbqt"};
  std::filesystem::path fld_path   = std::filesystem::path{"maps.fld"};

  po::options_description arguments_description("Available options");
  arguments_description.add_options()("help", "print this help message");
  arguments_description.add_options()("pdbqt",
                                      po::value(&pdbqt_path)->default_value(pdbqt_path),
                                      "Path to the protein file (in PDBQT)");
  arguments_description.add_options()("autogrid",
                                      po::value(&fld_path)->default_value(fld_path),
                                      "Path to the .fld file");
  // parse them
  po::options_description all("Allowed Options");
  all.add(arguments_description);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);

  po::notify(vm);

  auto protein_ptr = std::make_shared<mudock::dynamic_molecule>();
  auto& protein    = *protein_ptr;

  parse(protein, pdbqt_path);

  mudock::apply_autodock_forcefield_pdbqt(protein, pdbqt_path);

  //mudock::apply_autodock_forcefield(protein);
  const auto grid_atom_maps    = generate_atom_grid_maps(protein);
  const auto electrostatic_map = generate_electrostatic_grid_map(protein);
  const auto desolvation_map   = generate_desolvation_grid_map(protein);

  const auto desc = read_from_stream(std::ifstream(fld_path));
  std::stringstream desc_s{desc};

  std::array<mudock::ligand_map_types, mudock::num_ligand_map_types()> variables;
  int label       = 0;
  int label_eletr = 0;
  int label_desol = 0;
  std::string line;
  // Read the header (first few lines) for the grid information
  while (std::getline(desc_s, line)) {
    // Skip empty lines
    if (line.empty())
      continue;

    if (line.find(fld_tokens::LABEL_TOKEN) != std::string::npos) {
      size_t equal_pos = line.find('=');
      size_t dash_pos  = line.find('-');
      ++label;
      if (equal_pos != std::string::npos && dash_pos != std::string::npos && equal_pos < dash_pos) {
        const std::string result = line.substr(equal_pos + 1, dash_pos - equal_pos - 1);
        variables[label - 1]     = mudock::parse_map_symbol(result);
      }
      if (equal_pos != std::string::npos) {
        if (line.find(fld_tokens::ELETRO_TOKEN) != std::string::npos)
          label_eletr = label;
        else if (line.find(fld_tokens::DESOLV_TOKEN) != std::string::npos)
          label_desol = label;
      } else {
        throw std::runtime_error("Invalid label in fld");
      }
    }

    if (line.find(fld_tokens::VARIABLE_TOKEN) != std::string::npos &&
        line.find(fld_tokens::FILE_TOKEN) != std::string::npos) {
      int id;

      std::istringstream stream(line);
      std::string token, map_path;
      stream >> token >> id >> map_path >> token;
      const mudock::grid_map* reference_grid_map = &electrostatic_map;

      map_path = map_path.substr(std::strlen(fld_tokens::FILE_TOKEN));
      if (id == label_desol)
        reference_grid_map = &desolvation_map;
      else if (id != label_eletr && id >= 0 && id < mudock::num_ligand_map_types()) {
        reference_grid_map = &grid_atom_maps.get_atom_map(mudock::autodock_type_from_map(variables[id - 1]));
      } else if (id != label_eletr) {
        throw std::runtime_error("Unknown label/variables map in fld files");
      }
      mudock::grid_map autogrid_map = load_autogrid_map(map_path);
      // Compare eletrostatic map
      for (int k = 0; k < std::min(reference_grid_map->index.size_z(), autogrid_map.index.size_z()); ++k)
        for (int j = 0; j < std::min(reference_grid_map->index.size_y(), autogrid_map.index.size_y()); ++j)
          for (int i = 0; i < std::min(reference_grid_map->index.size_x(), autogrid_map.index.size_x()); ++i)
            if (std::abs(static_cast<float>(round3dp(reference_grid_map->at(i, j, k))) -
                         static_cast<float>(autogrid_map.at(i, j, k))) > float{0.01}) {
              mudock::error(std::format("Difference betweem maps {} at ({},{},{})", map_path, i, j, k));
              throw std::runtime_error("Error in Map");
            }
    }
  }
  mudock::info(std::format("Succesfully verified grid maps in {}", fld_path.string()));
  return EXIT_SUCCESS;
}
