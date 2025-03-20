#include <boost/program_options.hpp>
#include <cstdlib>
#include <filesystem>
#include <mudock/chem/ligand_maps.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/log.hpp>
#include <mudock/mudock.hpp>
#include <mudock/type_alias.hpp>
#include <stdexcept>
#include <string>
#include <tests/autogrid.hpp>

// utility function that reads the whole content of a stream
template<class stream_type>
inline auto read_from_stream(stream_type&& in) {
  assert(in.good());
  return std::string{std::istreambuf_iterator<std::string::value_type>{in},
                     std::istreambuf_iterator<std::string::value_type>{}};
}

int main(int argc, char* argv[]) {
  namespace po                       = boost::program_options;
  std::filesystem::path protein_path = std::filesystem::path{"protein.pdb"};
  std::filesystem::path fld_path     = std::filesystem::path{"maps.fld"};

  po::options_description arguments_description("Available options");
  arguments_description.add_options()("help", "print this help message");
  arguments_description.add_options()("protein",
                                      po::value(&protein_path)->default_value(protein_path),
                                      "Path to the protein file (in PDB)");
  arguments_description.add_options()("autogrid",
                                      po::value(&fld_path)->default_value(fld_path),
                                      "Path to the .fld file");
  // parse them
  po::options_description all("Allowed Options");
  all.add(arguments_description);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);

  po::notify(vm);

  mudock::info("Reading and parsing protein ", protein_path, " ...");
  auto protein_ptr = std::make_shared<mudock::dynamic_molecule>();
  auto& protein    = *protein_ptr;
  // auto pdb                       = mudock::pdb{};
  // const auto protein_description = read_from_stream(std::ifstream(protein_path));
  // pdb.parse(protein, protein_description);

  auto pdbqt                     = mudock::pdbqt{};
  const auto protein_description = read_from_stream(std::ifstream(protein_path));
  pdbqt.parse(protein, protein_description);

  mudock::apply_autodock_forcefield(protein);
  // auto grid_atom_maps    = std::make_shared<const mudock::grid_atom_mapper>(generate_atom_grid_maps(protein));
  auto electrostatic_map = std::make_shared<const mudock::grid_map>(generate_electrostatic_grid_map(protein));
  auto desolvation_map   = std::make_shared<const mudock::grid_map>(generate_desolvation_grid_map(protein));

  std::ifstream file(fld_path);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open file: " + fld_path.string());
  }
  std::string line;

  std::array<int, mudock::num_ligand_map_types()> variables;
  int label                = 0;
  int label_eletr          = 0;
  int label_desol          = 0;
  const std::string prefix = "file=";
  // Read the header (first few lines) for the grid information
  while (std::getline(file, line)) {
    // Skip empty lines
    if (line.empty())
      continue;

    if (line.find("label") != std::string::npos) {
      size_t equal_pos = line.find('=');
      size_t dash_pos  = line.find('-');
      ++label;
      if (equal_pos != std::string::npos && dash_pos != std::string::npos && equal_pos < dash_pos) {
        const std::string result = line.substr(equal_pos + 1, dash_pos - equal_pos - 1);
        variables[static_cast<int>(mudock::parse_map_symbol(result))] = label;
      }
      if (equal_pos != std::string::npos) {
        if (line.find("Electrostatics") != std::string::npos)
          label_eletr = label;
        else if (line.find("Desolvation") != std::string::npos)
          label_desol = label;

      } else {
        throw std::runtime_error("Invalid label in fld");
      }
    }

    if (line.find("variable") != std::string::npos && line.find(prefix) != std::string::npos) {
      int id;
      std::istringstream stream(line);
      std::string token, map_path;
      stream >> token >> id >> map_path >> token;
      if (id == label_eletr) {
        if (map_path.find(prefix) == 0) {
          map_path = map_path.substr(prefix.length());
        } else {
          throw std::runtime_error("Wrong eletrostatic map file");
        }

        GridMap eletro_map = loadGridMap(map_path);
        // Compare eletrostatic map
        for (int k = 0; k < eletro_map.z_size; ++k)
          for (int j = 0; j < eletro_map.y_size; ++j)
            for (int i = 0; i < eletro_map.z_size; ++i)
              if (std::abs(electrostatic_map->at(i, j, k) - eletro_map.get(i, j, k)) >
                  mudock::fp_type{0.005}) {
                mudock::error(std::format("Difference betweem eletrostatic map at ({},{},{}", i, j, k));
                throw std::runtime_error("Wrong Eletrostatic Map");
              }
      }

      if (id == label_desol) {
        if (map_path.find(prefix) == 0) {
          map_path = map_path.substr(prefix.length());
        } else {
          throw std::runtime_error("Wrong desolvation map file");
        }

        GridMap desolv_map = loadGridMap(map_path);
        // Compare eletrostatic map
        for (int k = 0; k < desolv_map.z_size; ++k)
          for (int j = 0; j < desolv_map.y_size; ++j)
            for (int i = 0; i < desolv_map.z_size; ++i)
              if (std::abs(desolvation_map->at(i, j, k) - desolv_map.get(i, j, k)) > mudock::fp_type{0.005}) {
                mudock::error(std::format("Difference betweem desolvation map at ({},{},{}", i, j, k));
                throw std::runtime_error("Wrong");
              }
      }
    }
  }

  return EXIT_SUCCESS;
}
