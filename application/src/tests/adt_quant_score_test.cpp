#include "autogrid.hpp"

#include <boost/program_options.hpp>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/log.hpp>
#include <mudock/mudock.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <stdexcept>
#include <string>
#include <mudock/chem/autodock_quant_protein.hpp>
#include <mudock/compute/adt_quant_score_kernel.hpp>

template<class T>
inline T round3dp(const T x) {
  return ((std::floor((x) * 1000.0 + 0.5)) / 1000.0);
}

std::vector<mudock::fp_type> generate_test_quantized_maps(const mudock::autodock_grid& adt_grid) {
    const auto& thresh = mudock::autodock_quant_protein::thresholds;
    int num_bins = thresh.size() + 1;
    size_t map_flat_size = adt_grid.get_size_xyz();
    
    std::vector<mudock::fp_type> quant_maps(num_bins * map_flat_size, 0.0);
    const mudock::fp_type* original_maps = adt_grid.get_maps_pointer();
    
    const mudock::fp_type* electro_map = original_maps + map_flat_size * static_cast<int>(mudock::autodock_grid_type::ELEC);
    const mudock::fp_type* desolv_map  = original_maps + map_flat_size * static_cast<int>(mudock::autodock_grid_type::DESOLV);

    for (int b = 0; b < num_bins; ++b) {
        mudock::fp_type charge_val;
        if (b == 0) charge_val = thresh[0] - 0.1f;
        else if (static_cast<size_t>(b) >= thresh.size()) charge_val = thresh.back() + 0.1f;
        else charge_val = (thresh[b] + thresh[b-1]) / 2;

        for (size_t i = 0; i < map_flat_size; ++i) {
            quant_maps[b * map_flat_size + i] = electro_map[i] * charge_val + desolv_map[i] * std::fabs(charge_val);
        }
    }
    return quant_maps;
}
// ----------------------------

int main(int argc, char *argv[]) {
  namespace po                   = boost::program_options;
  std::filesystem::path dpf_path = std::filesystem::path{"rec.dpf"};

  po::options_description arguments_description("Available options");
  arguments_description.add_options()("help", "print this help message");
  arguments_description.add_options()("dpf",
                                      po::value(&dpf_path)->default_value(dpf_path),
                                      "Path to the autodock DPF file");

  po::options_description all("Allowed Options");
  all.add(arguments_description);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
  po::notify(vm);

  mudock::autodock_grid adt_grid = load_autogrid_map_dpf(dpf_path);

  const auto ligand_path = get_ligand_path(dpf_path);
  mudock::static_molecule ligand =
      mudock::parser<mudock::static_molecule>(ligand_path, &mudock::pdbqt_rotate_check);

  auto f =
      std::function<void(mudock::autodock_static_layer &)>{[ligand_path](mudock::autodock_static_layer &l) {
        mudock::apply_autodock_forcefield_pdbqt(l, ligand_path);
      }};
  mudock::autodock_ligand adt_ligand{ligand, f};
  
  const auto adt_score       = load_autodock_score(dpf_path);
  const auto adt_error_score = load_autodock_error_score(dpf_path);

  const auto num_atoms    = ligand.num_atoms();
  const auto num_rotamers = ligand.num_rotamers();
  adt_ligand.update_offsets(static_cast<int>(adt_grid.get_map_flat_size()));

  mudock::info("Computing quantized FSR maps for custom kernel...");
  
  std::vector<mudock::fp_type> my_quant_maps = generate_test_quantized_maps(adt_grid);
  
  const auto& thresholds = mudock::autodock_quant_protein::thresholds;
  std::vector<int> atom_bins_b = mudock::build_atom_to_bin_map(ligand, thresholds);

  mudock::info("Computing energy ...");
  std::vector<int> num_atoms_b{num_atoms};
  std::vector<int> num_rotamers_b{static_cast<int>(num_rotamers)};
  std::vector<int> num_nonbonds_b{0, static_cast<int>(adt_ligand.non_bond_size())};
  std::vector<mudock::fp_type> scores_b{0};
  
  auto q = std::make_shared<mudock::queue_cpp>(0, mudock::device_type::CPU);
  mudock::adt_quant_score_kernel<mudock::queue_cpp> quant_kernel{ 
                                                         1,
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
                                                         my_quant_maps.data(),
                                                         atom_bins_b.data(), 
                                                         adt_grid.get_min_p(),
                                                         adt_grid.get_max_p(),
                                                         adt_grid.get_center_p(),
                                                         static_cast<int>(adt_grid.get_size_x()),
                                                         static_cast<int>(adt_grid.get_size_xy()),
                                                         static_cast<int>(adt_grid.get_size_xyz()),
                                                         scores_b.data(),
                                                         q};

  quant_kernel();
  const mudock::fp_type energy = scores_b[0];
  
  // Tollerance high because we are comparing a quantized score with a non-quantized one, so we expect some differences.
  mudock::info(std::format("Baseline Score: {}", adt_score - adt_error_score));
  mudock::info(std::format("Quantized Score: {}", energy));

 
  if (std::abs(energy - (adt_score - adt_error_score)) > mudock::fp_type{1.0}) { 
      mudock::error(std::format("Difference between scores of {} is too high ({} vs {})",
                                dpf_path.string(),
                                adt_score - adt_error_score,
                                energy));
      throw std::runtime_error("Error in score");
  }

  mudock::info(std::format("Successfully verified the QUANTIZED score in {}", dpf_path.string()));
  return EXIT_SUCCESS;
}