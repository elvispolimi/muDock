#include <boost/program_options.hpp>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>
#include <mudock/format.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/log.hpp>
#include <mudock/mudock.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <stdexcept>
#include <string>
#include <tests/autogrid.hpp>
#include <utility>

#ifndef MUDOCK_CTEST_SKIP_RETURN_CODE
#  error "MUDOCK_CTEST_SKIP_RETURN_CODE must be defined by CMake."
#endif

namespace {
// CTest treats this exit code as "skipped" when the grid test is not applicable.
constexpr int ctest_skip_return_code = MUDOCK_CTEST_SKIP_RETURN_CODE;
}

template<class T>
inline T round3dp(const T x) {
  return ((std::floor((x) *T{1000.0} + T{0.5})) / T{1000.0});
}

int main(int argc, char* argv[]) {
  if constexpr (std::is_same<mudock::fp_type, double>::value) {
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

    mudock::dynamic_molecule protein = mudock::parser<mudock::dynamic_molecule>(pdbqt_path);
    auto f =
        std::function<void(mudock::autodock_dynamic_layer&)>{[pdbqt_path](mudock::autodock_dynamic_layer& l) {
          mudock::apply_autodock_forcefield_pdbqt(l, pdbqt_path);
        }};
    mudock::autodock_protein adt_protein{protein, f};

    mudock::autodock_grid protein_autogrid = load_autogrid_map_fld(fld_path);

    for (int map_index = 0; map_index < mudock::num_autodock_grids(); ++map_index) {
      const auto map_type           = static_cast<mudock::autodock_grid_type>(map_index);
      const auto reference_grid_map = adt_protein.get_atom_map(map_type);
      const auto autogrid_map       = protein_autogrid.get_atom_map(map_type);

      for (size_t k = 0; k < std::min(reference_grid_map.z(), autogrid_map.z()); ++k)
        for (size_t j = 0; j < std::min(reference_grid_map.y(), autogrid_map.y()); ++j)
          for (size_t i = 0; i < std::min(reference_grid_map.x(), autogrid_map.x()); ++i) {
            const auto reference_round = static_cast<float>(round3dp(reference_grid_map.get(i, j, k)));
            const auto autogrid_round  = static_cast<float>(autogrid_map.get(i, j, k));
            const auto max_absolute = std::max(std::fabs(reference_round), std::fabs(autogrid_round)) / 100;
            const auto delta        = std::clamp(max_absolute, float{0.01}, float{1});
            if (std::fabs(reference_round - autogrid_round) > delta) {
              mudock::error(std::format(
                  "Difference betweem maps {} at ({},{},{}): muDock {} autogrid {} with an error threshold of {}",
                  mudock::get_description(map_type).name,
                  i,
                  j,
                  k,
                  reference_round,
                  autogrid_round,
                  delta));
              throw std::runtime_error("Error in Map");
            }
          }
    }
    mudock::info(std::format("Succesfully verified grid maps in {}", fld_path.string()));
  } else {
    mudock::info("Grid test requires mudock::fp_type to be double");
    return ctest_skip_return_code;
  }
  return EXIT_SUCCESS;
}
