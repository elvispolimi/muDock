

// #include <boost/program_options.hpp>
// #include <cstring>
// #include <filesystem>
// #include <fstream>
// #include <functional>
// #include <iostream>
// #include <memory>
// #include <string>
// #include <vector>

// #include <mudock/batch.hpp>
// #include <mudock/chem/autodock_grid_types.hpp>
// #include <mudock/compute/adadelta.hpp>
// #include <mudock/compute/adt_score.hpp>
// #include <mudock/compute/geometric_transform.hpp>
// #include <mudock/compute/scratchpad.hpp>
// #include <mudock/cpp_implementation/queue_cpp.hpp>
// #include <mudock/format/reader.hpp>
// #include <mudock/format/supported_format.hpp>
// #include <mudock/log.hpp>
// #include <mudock/molecule.hpp>
// #include <mudock/mudock.hpp>
// #include <mudock/type_alias.hpp>

// struct local_search_options {
//   std::filesystem::path protein_path;
//   std::filesystem::path ligand_path;
//   std::string ls_method = "adadelta";
//   std::size_t population = 1;
//   std::size_t ls_iterations = 300;
// };

// local_search_options parse_local_search_options(int argc, char* argv[]) {
//   namespace po = boost::program_options;
//   local_search_options options;

//   po::options_description description("Local search options");
//   description.add_options()("help", "Print this help message");
//   description.add_options()("protein",
//                              po::value(&options.protein_path)->required(),
//                              "Protein file path (PDBQT/PDB)");
//   description.add_options()("ligand",
//                              po::value(&options.ligand_path)->required(),
//                              "Ligand file path (MOL2/PDBQT/PDB)");
//   description.add_options()("ls_method",
//                              po::value(&options.ls_method)->default_value(options.ls_method),
//                              "Local search method (currently only 'adadelta')");
//   description.add_options()("population",
//                              po::value(&options.population)->default_value(options.population),
//                              "Number of individuals per ligand for local search");
//   description.add_options()("ls_iterations",
//                              po::value(&options.ls_iterations)->default_value(options.ls_iterations),
//                              "Number of local search iterations");

//   po::variables_map vm;
//   po::store(po::command_line_parser(argc, argv).options(description).run(), vm);

//   if (vm.count("help") > 0) {
//     std::cout << "Standalone local-search executable for a single protein-ligand pair." << std::endl
//               << std::endl;
//     std::cout << "USAGE: " << argv[0]
//               << " --protein <protein.pdbqt> --ligand <ligand.mol2> [--ls_method adadelta]"
//               << std::endl
//               << std::endl;
//     std::cout << description << std::endl;
//     std::exit(EXIT_SUCCESS);
//   }

//   po::notify(vm);
//   return options;
// }

// template<typename queue_t>
// mudock::fp_type read_score_from_scratch(std::shared_ptr<mudock::scratchpad<queue_t>> scratch) {
//   auto &scores = (*scratch).template get<mudock::buffer_data_type::SCORES>();
//   scores.copy_device2host();
//   (*scratch).get_queue()->synchronize();
//   return scores()[0];
// }

// int main(int argc, char* argv[]) {
//   const auto options = parse_local_search_options(argc, argv);
//   mudock::info("Running local search on protein=", options.protein_path, " ligand=", options.ligand_path);
//   mudock::info("Local search method=", options.ls_method, " population=", options.population,
//                 " iterations=", options.ls_iterations);

//   const auto protein = std::make_shared<mudock::dynamic_molecule>(
//       mudock::parser<mudock::dynamic_molecule>(options.protein_path));

//   auto ligand = std::make_unique<mudock::static_molecule>(
//       mudock::parser<mudock::static_molecule>(options.ligand_path));

//   mudock::batch<mudock::static_molecule> batch;
//   batch.num_ligands        = 1;
//   batch.batch_max_atoms    = ligand->num_atoms();
//   batch.batch_max_rotamers = ligand->num_rotamers();
//   batch.molecules[0]       = std::move(ligand);

//   mudock::knobs knobs;
//   knobs.population_number = std::max<std::size_t>(1, options.population);
//   knobs.ls_iterations     = options.ls_iterations;

//   using queue_t = mudock::queue_cpp;
//   auto scratch = std::make_shared<mudock::scratchpad<queue_t>>(knobs, 0, mudock::device_type::CPU);
//   auto device_scratch = std::make_shared<mudock::scratchpad<queue_t>>(knobs, 0, mudock::device_type::CPU);

//   auto score_stage = std::make_shared<mudock::adt_score<queue_t>>(scratch, device_scratch, *protein);
//   mudock::geometric<queue_t> geom_transform(scratch, *protein);

//   std::unique_ptr<mudock::local_search<queue_t, mudock::adt_score>> local_search_stage;
//   if (options.ls_method == "adadelta") {
//     local_search_stage = std::make_unique<mudock::adadelta<queue_t, mudock::adt_score>>(scratch, score_stage);
//   } else {
//     mudock::error("Unsupported local search method: ", options.ls_method);
//   }

//   auto &chromosomes = (*scratch).template get<mudock::buffer_data_type::CHROMOSOMES>();
//   chromosomes.alloc(static_cast<std::size_t>(knobs.population_number) * batch.num_ligands);
//   if (chromosomes.num_elements() > 0) {
//     std::memset(chromosomes(), 0, chromosomes.num_elements() * sizeof(mudock::chromosome));
//   }
//   chromosomes.set_valid();

//   geom_transform.prepare(batch);
//   score_stage->prepare(batch);
//   local_search_stage->prepare(batch);

//   score_stage->operator()();
//   const auto initial_score = read_score_from_scratch(scratch);
//   mudock::info("Initial score=", initial_score);

//   if (auto* adadelta_ptr = dynamic_cast<mudock::adadelta<queue_t, mudock::adt_score> *>(
//           local_search_stage.get())) {
//     adadelta_ptr->set_coordinate_update([&]() { geom_transform(); });
//   }
  
//   // non ho capito perché geom transform va messo così (local stage operator e poi geom transform), altrimenti:
//   // - se lo metto prima, le ls_iter non fanno cambiare il risultato
//   // - se lo tolgo, spara sempre numeroni come se non si inizializzasse 
//   local_search_stage->operator()();
//   geom_transform();
//   score_stage->operator()();
//   score_stage->teardown(batch);
//   geom_transform.teardown(batch);

//   std::cout << "Ligand " << batch.molecules[0]->properties.get(mudock::property_type::NAME)
//             << " final score " << batch.molecules[0]->properties.get(mudock::property_type::SCORE)
//             << std::endl;

//   mudock::info("Local search completed.");
//   return EXIT_SUCCESS;
// }





#include "command_line_args.hpp"

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <fstream>
#include <memory>
#include <mutex>
#include <mudock/compute/manager.hpp>
#include <mudock/compute/pipeline.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>
#include <thread>

int main(int argc, char* argv[]) {
  const auto args = parse_command_line_arguments(argc, argv);

  mudock::info("Reading and parsing protein ", args.protein_path, " ...");
  auto protein =
      std::make_shared<mudock::dynamic_molecule>(mudock::parser<mudock::dynamic_molecule>(args.protein_path));

  mudock::info("Reading ligand ", args.ligand_path, " ...");
  const auto in_format = mudock::parse_supported_format(args.ligand_path);
  auto input_queue     = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
  constexpr_switch<0, mudock::get_num_supported_format(), 1>(
      [&](const auto format_index) {
        const auto format = static_cast<mudock::supported_format>(format_index());
        auto input_text   = read_from_stream(std::ifstream(args.ligand_path));
        mudock::splitter<mudock::type_of_format<static_cast<mudock::supported_format>(format_index())>> split;
        auto ligands_description = split(std::move(input_text));
        ligands_description.emplace_back(split.flush());
        input_queue->initialize(ligands_description.size());

        mudock::info("Parsing ", ligands_description.size(), " ligand(s) ...");
        if constexpr (format == mudock::supported_format::ADTMOL2) {
#ifdef _OPENMP
#pragma omp parallel for shared(input_queue)
#endif
          for (const auto& description: ligands_description) {
            auto ligand = std::make_unique<mudock::static_molecule>(
                mudock::parser<mudock::supported_format::ADTMOL2, mudock::static_molecule>(description));
            input_queue->enqueue(ligand);
          }
        } else {
          for (const auto& description: ligands_description) {
            try {
              auto ligand = std::make_unique<mudock::static_molecule>(
                  mudock::parser<format, mudock::static_molecule>(description));
              input_queue->enqueue(ligand);
            } catch (...) {}
          }
        }
      },
      in_format);

  input_queue->send_terminate_signal(); // signal that no more ligand will be enqueued in the input queue    

  mudock::info("Running local search...");
  // mudock::info("Scores per ligand (population): ", args.knobs.population_number);

  mudock::local_search_pipeline pipe{protein};
  auto output_queue = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
  output_queue->initialize(input_queue->size());

  const auto start = std::chrono::high_resolution_clock::now();
  std::atomic<std::size_t> dropped_by_timeout{0};
  std::atomic<bool> timeout_triggered{false};
  {
    std::mutex observer_mutex;
    std::condition_variable observer_cv;
    bool observer_stop = false;
    std::thread observer_thread;
    if (args.observer && *args.observer > 0.0) {
      mudock::info("Observer enabled with period: ", *args.observer, " s");
      observer_thread = std::thread([&]() {
        std::size_t prev_processed = output_queue->size();
        auto prev_time             = std::chrono::high_resolution_clock::now();
        while (true) {
          std::unique_lock<std::mutex> lock(observer_mutex);
          const bool stop = observer_cv.wait_for(lock,
                                                 std::chrono::duration<double>(*args.observer),
                                                 [&]() { return observer_stop; });
          if (stop) {
            break;
          }
          lock.unlock();

          const auto now                  = std::chrono::high_resolution_clock::now();
          const std::size_t now_processed = output_queue->size();
          const std::size_t in_backlog    = input_queue->size();

          const std::chrono::duration<double> dt = now - prev_time;
          const std::size_t delta_processed       = now_processed - prev_processed;
          const double inst_throughput =
              dt.count() > 0.0 ? static_cast<double>(delta_processed) / dt.count() : 0.0;
          const std::chrono::duration<double> total = now - start;
          const double avg_throughput =
              total.count() > 0.0 ? static_cast<double>(now_processed) / total.count() : 0.0;

          mudock::info("Observer: processed=",
                       now_processed,
                       ", input_backlog=",
                       in_backlog,
                       ", inst_throughput=",
                       inst_throughput,
                       " ligands/s, avg_throughput=",
                       avg_throughput,
                       " ligands/s");

          prev_processed = now_processed;
          prev_time      = now;
        }
      });
    }

    std::mutex timer_mutex;
    std::condition_variable timer_cv;
    bool timer_cancelled = false;
    std::thread timer_thread;
    if (args.time_limit_sec && *args.time_limit_sec > 0.0) {
      mudock::info("Time limit enabled: ", *args.time_limit_sec, " s");
      timer_thread = std::thread([&]() {
        std::unique_lock<std::mutex> lock(timer_mutex);
        const bool cancelled = timer_cv.wait_for(lock,
                                                 std::chrono::duration<double>(*args.time_limit_sec),
                                                 [&]() { return timer_cancelled; });
        if (cancelled) {
          return;
        }
        lock.unlock();
        input_queue->send_terminate_signal();
        dropped_by_timeout.store(input_queue->clear(), std::memory_order_relaxed);
        timeout_triggered.store(true, std::memory_order_relaxed);
        mudock::info("Time limit reached: discarded ",
                     dropped_by_timeout.load(std::memory_order_relaxed),
                     " pending ligand(s) from input queue.");
      });
    }

    {
      auto threadpool = mudock::threadpool();
      mudock::manager(args.device_confs, threadpool, args.knobs, input_queue, output_queue, pipe);
    } // threadpool destructor waits for workers; computation is complete here

    output_queue->send_terminate_signal(); // signal that no more ligand will be enqueued in the output queue
    if (observer_thread.joinable()) {
      {
        std::lock_guard<std::mutex> lock(observer_mutex);
        observer_stop = true;
      }
      observer_cv.notify_one();
      observer_thread.join();
    }

    if (timer_thread.joinable()) {
      {
        std::lock_guard<std::mutex> lock(timer_mutex);
        timer_cancelled = true;
      }
      timer_cv.notify_one();
      timer_thread.join();
    }
  }
  

  const auto end                              = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> elapsed = end - start;

  // std::size_t processed = 0;
  // for (auto ligand = output_queue->dequeue(); ligand; ligand = output_queue->dequeue()) { ++processed; }
  auto ligand = output_queue->dequeue();
  std::cout << ligand->properties.get(mudock::property_type::NAME) << " "
            << ligand->properties.get(mudock::property_type::SCORE) << std::endl;

  if (timeout_triggered.load(std::memory_order_relaxed)) {
    mudock::info("Dropped ligands due to timeout: ", dropped_by_timeout.load(std::memory_order_relaxed));
  }
  // mudock::info("Elapsed time: ", elapsed.count(), " s");
  // if (elapsed.count() > 0.0) {
  //   mudock::info("Throughput: ", static_cast<double>(processed) / elapsed.count(), " ligands/s");
  // }

  mudock::info("Local search completed.");
  return EXIT_SUCCESS;
}
