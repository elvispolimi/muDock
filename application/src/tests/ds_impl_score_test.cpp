#include <boost/program_options.hpp>
#include <command_line_args.hpp>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/drugscore/ds_score.hpp>
#include <mudock/compute/manager.hpp>
#include <mudock/compute/pipeline.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/knobs.hpp>
#include <mudock/log.hpp>
#include <mudock/mudock.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <sstream>
#include <string>
#include <tests/autogrid.hpp>
#include <utility>

template<typename pipeline_t>
mudock::fp_type run_config(const mudock::knobs& knobs,
                           std::shared_ptr<mudock::safe_queue<mudock::static_molecule>>& input_queue,
                           std::shared_ptr<mudock::safe_queue<mudock::static_molecule>>& output_queue,
                           std::shared_ptr<mudock::static_molecule>& ligand,
                           pipeline_t& pipe,
                           const std::string_view& use_config
                          ) {
  input_queue->clear_terminate_signal(); // clear the terminate signal to be able to enqueue new ligands
  output_queue->clear_terminate_signal();

  auto ligand_in    = std::make_unique<mudock::static_molecule>(*ligand);
  input_queue->enqueue(ligand_in);
  input_queue->send_terminate_signal(); // signal that no more ligand will be enqueued in the input queue
  {
    auto threadpool = mudock::threadpool();
    mudock::manager({std::string{use_config}}, threadpool, knobs, input_queue, output_queue, pipe);
  }
  output_queue->send_terminate_signal(); // signal that no more ligand will be enqueued in the output queue
  auto ligand_out = output_queue->dequeue();
  std::stringstream ss(ligand_out->properties.get(mudock::property_type::SCORE));
  mudock::fp_type score;
  ss >> score;

  return score;
}


int main(int argc, char* argv[]) {
  namespace po                          = boost::program_options;
  std::filesystem::path protein_path    = std::filesystem::path{"protein.pdb"};
  std::filesystem::path ligand_path     = std::filesystem::path{"ligand.mol2"};
  const auto knobs                      = mudock::knobs{1, 1, 0, 0, 123};
  std::vector<std::string> device_confs = {std::string{use_cpu_conf}};

  po::options_description arguments_description("Available options");
  arguments_description.add_options()("help", "print this help message");
  arguments_description.add_options()("protein,p",
                                      po::value(&protein_path)->default_value(protein_path),
                                      "Path to the protein file");
  arguments_description.add_options()("ligand,l",
                                      po::value(&ligand_path)->default_value(ligand_path),
                                      "Path to the ligand file");
  arguments_description.add_options()(
      "use",
      po::value<std::vector<std::string>>(&device_confs)->multitoken()->composing(),
      "Map each implementation to the device");

  // parse them
  po::options_description all("Allowed Options");
  all.add(arguments_description);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);

  po::notify(vm);

  mudock::info("Reading and parsing protein ", protein_path, " ...");
  auto protein =
      std::make_shared<mudock::dynamic_molecule>(mudock::parser<mudock::dynamic_molecule>(protein_path));

  mudock::info("Reading and parsing ligand ", ligand_path, " ...");
  auto ligand =
      std::make_shared<mudock::static_molecule>(mudock::parser<mudock::static_molecule>(ligand_path));

  mudock::ds_score_pipeline pipe{protein};
  mudock::ds_precomp_score_pipeline pipe_precomp{protein};

  mudock::info("Generating score reference ...");
  auto output_queue = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
  auto input_queue  = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();

  mudock::fp_type reference_score = run_config(knobs, input_queue, output_queue, ligand, pipe, use_cpu_conf);
  mudock::fp_type reference_score_precomp = run_config(knobs, input_queue, output_queue, ligand, pipe_precomp, use_cpu_conf);
  

  for (auto& conf: device_confs) {
    mudock::info("Comparing reference with ", conf, " ...");
    
    mudock::fp_type score = run_config(knobs, input_queue, output_queue, ligand, pipe, conf);
    mudock::fp_type score_precomp = run_config(knobs, input_queue, output_queue, ligand, pipe_precomp, conf);

    const auto diff = std::abs(score - reference_score);
    const auto diff_precomp = std::abs(score_precomp - reference_score_precomp);
    
    const auto max_absolute = std::max(std::fabs(score), std::fabs(reference_score)) / 1000;
    const auto max_absolute_precomp = std::max(std::fabs(score_precomp), std::fabs(reference_score_precomp)) / 1000;

    const auto error        = std::clamp(max_absolute, mudock::fp_type{0.001f}, mudock::fp_type{5});
    const auto error_precomp = std::clamp(max_absolute_precomp, mudock::fp_type{0.001f}, mudock::fp_type{5});

    if (diff > error) {
      mudock::error(std::format(
          "Difference betweem scores of {} on {} ( CPU {} vs {} {} with an error threshold of {})",
          ligand_path.string(),
          protein_path.string(),
          reference_score,
          conf,
          score,
          error));

      return EXIT_FAILURE;
    }

    if (diff_precomp > error_precomp) {
      mudock::error(std::format(
          "Difference betweem scores with precompute of {} on {} ( CPU {} vs {} {} with an error threshold of {})",
          ligand_path.string(),
          protein_path.string(),
          reference_score_precomp,
          conf,
          score_precomp,
          error_precomp));

      return EXIT_FAILURE;
    }
  }

  mudock::info(std::format("Succesfully verified the score {} and score precompute {} of {} on {} and ",
                           reference_score,
                           reference_score_precomp,
                           ligand_path.string(),
                           protein_path.string()));
  return EXIT_SUCCESS;
}
