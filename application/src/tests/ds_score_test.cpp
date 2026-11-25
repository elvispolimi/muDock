#include <boost/program_options.hpp>
#include <chrono>
#include <command_line_args.hpp>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/likwid_utils.hpp>
#include <mudock/log.hpp>
#include <mudock/mudock.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <regex>
#include <optional>
#include <sstream>
#include <string>

using namespace mudock;

std::optional<fp_type> readDrugScore(const std::filesystem::path& path){
    std::ifstream file(path);
    if (!file)
        return std::nullopt;

    const std::regex pattern(
        R"(^##\s+DRUGSCORE\s+([-+]?(?:\d*\.?\d+)(?:[eE][-+]?\d+)?)\s*$)"
    );

    std::string line;
    std::smatch match;

    while (std::getline(file, line))
    {
        if (std::regex_match(line, match, pattern))
        {
            return std::stod(match[1].str());
        }
    }

    return std::nullopt;
}

int main(int argc, char* argv[]) {
  // TODO: Reformat on new call method of score function
  namespace po                       = boost::program_options;
  std::filesystem::path protein_path = std::filesystem::path{"protein.pdb"};
  std::filesystem::path ligand_path  = std::filesystem::path{"ligand.mol2"};
  std::filesystem::path ref_score_path = std::filesystem::path("ligand_score.ds");

  po::options_description arguments_description("Available options");
  arguments_description.add_options()("help", "print this help message");
  arguments_description.add_options()("protein,p",
                                      po::value(&protein_path)->default_value(protein_path),
                                      "Path to the protein file");
  arguments_description.add_options()("ligand,l",
                                      po::value(&ligand_path)->default_value(ligand_path),
                                      "Path to the ligand file");
  arguments_description.add_options()("score,s",
                                      po::value(&ref_score_path)->default_value(ref_score_path),
                                      "Path to the ligand reference score file");

  // parse them
  po::options_description all("Allowed Options");
  all.add(arguments_description);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);

  po::notify(vm);

  mudock::info("Reading and parsing protein ", protein_path, " ...");
  auto protein = mudock::parser<mudock::dynamic_molecule>(protein_path, &mudock::pdbqt_rotate_check);

  mudock::drugscore_protein ds_prot(protein);
  mudock::drugscore_precomp_protein ds_prot_precomp(protein);


  mudock::info("Reading and parsing ligand ", ligand_path, " ...");
  const auto ligand_format = mudock::parse_supported_format(std::filesystem::path{ligand_path});
  const auto rotate_check  = ligand_format == mudock::supported_format::PDBQT ? &mudock::pdbqt_rotate_check
                                                                               : &mudock::ob_rotate_check;
  
  mudock::static_molecule ligand = mudock::parser<mudock::static_molecule>(ligand_path, rotate_check);   
  mudock::drugscore_ligand ds_lig{ligand};
  auto ref_score = readDrugScore(ref_score_path);
  if (!ref_score){
    mudock::info(std::format("Failed test of {} on {}: reference score not found",
                            ligand_path.string(),
                            protein_path.string())); 
                          
    return EXIT_FAILURE;
  }

  const fp_type* __restrict__ ligand_x = ds_lig.get_base_molecule().x();
  const fp_type* __restrict__ ligand_y = ds_lig.get_base_molecule().y();
  const fp_type* __restrict__ ligand_z = ds_lig.get_base_molecule().z();
  const int num_atoms                  = ds_lig.get_base_molecule().num_atoms();
  std::vector<int> num_atoms_b{num_atoms};
  const ds_type* __restrict__ ligand_ds_types = ds_lig.drugscore_type();

  const fp_type* __restrict__ protein_x        = ds_prot.get_base_molecule().x();
  const fp_type* __restrict__ protein_y        = ds_prot.get_base_molecule().y();
  const fp_type* __restrict__ protein_z        = ds_prot.get_base_molecule().z();
  const int num_protein_atoms                  = ds_prot.get_base_molecule().num_atoms();
  const ds_type* __restrict__ protein_ds_types = ds_prot.drugscore_type();

  const fp_type* __restrict__ ds_scores = DRUGSCORE_MODEL_FLAT();

  const fp_type* __restrict__ prot_min          = ds_prot_precomp.get_min_p();
  const fp_type* __restrict__ prot_max          = ds_prot_precomp.get_max_p();
  const int prot_index_x                        = static_cast<int>(ds_prot_precomp.get_size_x());
  const int prot_index_xy                       = static_cast<int>(ds_prot_precomp.get_size_xy());
  const int prot_index_xyz                      = static_cast<int>(ds_prot_precomp.get_size_xyz());
  const fp_type* __restrict__ prot_ds_grid_maps = ds_prot_precomp.get_maps_pointer();

  std::vector<mudock::fp_type> scores{0};
  auto q = std::make_shared<mudock::queue_cpp>(0, mudock::device_type::CPU);

  mudock::ds_score_kernel<mudock::queue_cpp> ds_kernel{num_atoms,
                                                        1,
                                                        1,
                                                        ligand_x,
                                                        ligand_y,
                                                        ligand_z,
                                                        ligand_ds_types,
                                                        num_atoms_b.data(),
                                                        protein_x,
                                                        protein_y,
                                                        protein_z,
                                                        protein_ds_types,
                                                        num_protein_atoms,
                                                        ds_scores,
                                                        scores.data(),
                                                        q};

  mudock::ds_score_precomp_kernel<mudock::queue_cpp> ds_kernel_precomp(num_atoms,
                                                                1,
                                                                1,
                                                                ligand_x,
                                                                ligand_y,
                                                                ligand_z,
                                                                ligand_ds_types,
                                                                num_atoms_b.data(),
                                                                prot_min,
                                                                prot_max,
                                                                prot_index_x,
                                                                prot_index_xy,
                                                                prot_index_xyz,
                                                                prot_ds_grid_maps,
                                                                scores.data(),
                                                                q);

  mudock::info("Scoring ...");
  MUDOCK_MARKER_INIT
  MUDOCK_CPP_MARKER_REGISTER("CPP ds_score_kernel");
  MUDOCK_CPP_MARKER_START("CPP ds_score_kernel");
  ds_kernel();
  MUDOCK_CPP_MARKER_STOP("CPP ds_score_kernel");

  const auto drugscore = scores[0];

  MUDOCK_CPP_MARKER_REGISTER("CPP ds_score_kernel");
  MUDOCK_CPP_MARKER_START("CPP ds_score_kernel");
  ds_kernel_precomp();
  MUDOCK_CPP_MARKER_STOP("CPP ds_score_kernel");
  MUDOCK_MARKER_CLOSE

  const auto drugscore_precomp = scores[0];

  
  fp_type threshold  = std::abs(*ref_score / 10); // static_cast<mudock::fp_type>(0.1)

  if (std::abs(drugscore - *ref_score) > threshold) {
    mudock::error(std::format("Difference betweem scores of {} ({} vs {})",
                              ligand_path.string(),
                              drugscore,
                              *ref_score));
    return EXIT_FAILURE;
  }
  if (std::abs(drugscore_precomp - *ref_score) > threshold) {
    mudock::error(std::format("Difference betweem scores precompute of {} ({} vs {})",
                              ligand_path.string(),
                              drugscore_precomp,
                              *ref_score));
    return EXIT_FAILURE;
  }

  mudock::info(std::format("Succesfully verified the score {} and score precompute {} of {} on {}",
                            drugscore,
                            drugscore_precomp,
                            ligand_path.string(),
                            protein_path.string()));

  return EXIT_SUCCESS;
}
