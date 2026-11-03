#include <mudock/compute/adadelta.hpp>
#include <mudock/cpp_implementation/adadelta_convergence_cpp.hpp>
#include <mudock/cpp_implementation/adadelta_cpp.hpp>
#include <mudock/cpp_implementation/geom_transform_cpp.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
  using namespace mudock;

  void require(const bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
  }

  chromosome pose(const fp_type value) {
    chromosome result;
    result.fill(value);
    return result;
  }

  struct kernel_fixture {
    std::vector<int> converged;
    std::vector<int> selected;
    std::vector<fp_type> scores;
    std::vector<fp_type> best_scores;
    std::vector<int> ages;
    std::vector<chromosome> current;
    std::vector<chromosome> best;
    std::shared_ptr<queue_cpp> q = std::make_shared<queue_cpp>(0, device_type::CPU);
    adadelta_convergence_kernel<queue_cpp> kernel;

    kernel_fixture(const int ligands = 1, const int population = 1)
        : converged(ligands, 0),
          selected(ligands * population, 1),
          scores(ligands * population, fp_type{10}),
          best_scores(ligands * population, std::numeric_limits<fp_type>::infinity()),
          ages(ligands * population, 0),
          current(ligands * population, pose(fp_type{1})),
          best(ligands * population),
          kernel(ligands, population, converged.data(), selected.data(), scores.data(),
                 best_scores.data(), ages.data(), current.data(), best.data(), q) {}
  };

  void warmup_and_best_pose() {
    kernel_fixture f;
    f.kernel(0, false);
    require(f.best[0] == pose(fp_type{1}), "Initial score must snapshot every chromosome gene");
    for (std::size_t iteration = 1; iteration < MIN_ITERATIONS_TO_WARM_UP; ++iteration) {
      f.current[0] = pose(static_cast<fp_type>(iteration + 1));
      f.kernel(iteration, false);
      require(f.selected[0] == 1, "Individual stopped before completing warmup");
    }
    f.current[0] = pose(fp_type{99});
    f.kernel(MIN_ITERATIONS_TO_WARM_UP, false);
    require(f.selected[0] == 0, "Stalled individual must stop exactly at the warmup boundary");
    require(f.current[0] == pose(fp_type{1}), "Stopping must restore the best complete chromosome");
  }

  void patience_and_improvement() {
    kernel_fixture f;
    f.kernel(MIN_ITERATIONS_TO_WARM_UP, false);
    for (int age = 1; age < MAX_BEST_SCORE_AGE; ++age) {
      f.kernel(MIN_ITERATIONS_TO_WARM_UP + static_cast<std::size_t>(age), false);
      require(f.selected[0] == 1, "Individual stopped before reaching patience threshold");
    }
    f.scores[0] = fp_type{9};
    f.current[0] = pose(fp_type{2});
    f.kernel(MIN_ITERATIONS_TO_WARM_UP + MAX_BEST_SCORE_AGE, false);
    require(f.ages[0] == 0 && f.selected[0] == 1, "Improvement must reset patience");
    f.current[0] = pose(fp_type{3});
    f.scores[0] = fp_type{11};
    for (int age = 1; age <= MAX_BEST_SCORE_AGE; ++age) {
      f.kernel(MIN_ITERATIONS_TO_WARM_UP + MAX_BEST_SCORE_AGE + static_cast<std::size_t>(age), false);
      require(f.selected[0] == (age < MAX_BEST_SCORE_AGE), "Incorrect patience threshold");
    }
    require(f.best_scores[0] == fp_type{9} && f.current[0] == pose(fp_type{2}),
            "Regression after improvement must restore the improved pose");
  }

  void individual_and_ligand_isolation() {
    kernel_fixture f(3, 2);
    f.selected[3] = 0;
    f.converged[2] = 7;
    f.kernel(0, false);
    f.ages[0] = f.ages[1] = f.ages[2] = MAX_BEST_SCORE_AGE - 1;
    f.scores[1] = fp_type{9};
    std::fill(f.current.begin(), f.current.end(), pose(fp_type{2}));
    f.kernel(MIN_ITERATIONS_TO_WARM_UP, false);
    require(f.selected == std::vector<int>({0, 1, 0, 0, 1, 1}),
            "Convergence must be tracked separately for each active individual");
    require(f.ages[1] == 0 && f.best_scores[1] == fp_type{9},
            "One individual's improvement must not affect its neighbors");
    require(f.ages[3] == 0 && std::isinf(f.best_scores[3]),
            "Unselected individual must remain untouched");
    require(f.ages[4] == 0 && f.ages[5] == 0 && std::isinf(f.best_scores[4]),
            "Already converged ligand must remain untouched");
    require(f.current[0] == pose(fp_type{1}) && f.current[1] == pose(fp_type{2}) &&
                f.current[2] == pose(fp_type{1}) && f.current[4] == pose(fp_type{2}),
            "Restoring a pose must not modify other individuals or ligands");
  }

  void final_score_and_nonfinite_scores() {
    kernel_fixture f;
    f.kernel(0, false);
    f.scores[0] = fp_type{5};
    f.current[0] = pose(fp_type{2});
    f.kernel(1, true);
    require(f.best_scores[0] == fp_type{5} && f.current[0] == pose(fp_type{2}),
            "The last AdaDelta update must be scored before final restoration");

    kernel_fixture worse;
    worse.kernel(0, false);
    worse.scores[0] = fp_type{20};
    worse.current[0] = pose(fp_type{2});
    worse.kernel(1, true);
    require(worse.current[0] == pose(fp_type{1}), "Finalization must restore an earlier better pose");

    kernel_fixture invalid;
    invalid.scores[0] = std::numeric_limits<fp_type>::quiet_NaN();
    invalid.kernel(0, false);
    invalid.scores[0] = -std::numeric_limits<fp_type>::infinity();
    invalid.current[0] = pose(fp_type{3});
    invalid.kernel(1, true);
    require(std::isinf(invalid.best_scores[0]) && invalid.best_scores[0] > 0,
            "Nonfinite scores must not establish a best pose");
    require(invalid.current[0] == pose(fp_type{3}),
            "An individual without a finite score must not restore an uninitialized pose");
  }

  void stage_reset_and_pingpong() {
    knobs conf;
    conf.population_number = 1;
    auto scratch = std::make_shared<scratchpad<queue_cpp>>(conf, 0, device_type::CPU);
    initialize_buffer<buffer_data_type::CONVERGED_LIGANDS>(scratch, 1, 0);
    initialize_buffer<buffer_data_type::LS_SELECTED>(scratch, 1, 1);
    initialize_buffer<buffer_data_type::SCORES>(scratch, 1, fp_type{10});
    initialize_buffer<buffer_data_type::CHROMOSOMES>(scratch, 1, pose(fp_type{1}));
    batch<static_molecule> batch;
    batch.num_ligands = 1;
    adadelta_convergence<queue_cpp> stage(scratch);
    stage.prepare(batch);
    auto& selected = scratch->get<buffer_data_type::LS_SELECTED>();
    auto& scores = scratch->get<buffer_data_type::SCORES>();
    auto& current = scratch->get<buffer_data_type::CHROMOSOMES>();
    auto& best_scores = scratch->get<buffer_data_type::BEST_SCORES>();
    auto& ages = scratch->get<buffer_data_type::BEST_SCORES_AGE>();
    const auto best_pointer = best_scores.dev_pointer();
    const auto age_pointer = ages.dev_pointer();
    const auto pose_pointer = scratch->get<buffer_data_type::ADADELTA_BEST_CHROMOSOMES>().dev_pointer();
    stage.reset();
    for (int iteration = 0; iteration < MIN_ITERATIONS_TO_WARM_UP; ++iteration) stage();
    require(selected()[0] == 1, "Stage warmup counter started before the first completed update");
    current()[0] = pose(fp_type{2});
    stage.finalize();
    require(selected()[0] == 0 && current()[0] == pose(fp_type{1}),
            "Finalization must count the last completed update and restore the best pose");

    selected()[0] = 1;
    scores()[0] = fp_type{100};
    std::vector<chromosome> alternate(1, pose(fp_type{4}));
    stage.set_chromosomes_buffer(alternate.data());
    stage.reset();
    require(best_scores.dev_pointer() == best_pointer && ages.dev_pointer() == age_pointer &&
                scratch->get<buffer_data_type::ADADELTA_BEST_CHROMOSOMES>().dev_pointer() == pose_pointer,
            "Generation reset must preserve pointers held by prepared kernels");
    require(ages()[0] == 0 && best_scores()[0] > scores()[0],
            "Generation reset must discard previous best score and patience");
    stage();
    require(selected()[0] == 1, "Generation reset must restart the warmup counter");
    alternate[0] = pose(fp_type{5});
    scores()[0] = fp_type{101};
    stage.finalize();
    require(alternate[0] == pose(fp_type{4}) && current()[0] == pose(fp_type{1}),
            "Convergence must follow the current genetic ping-pong chromosome buffer");
  }

  template<typename queue_t>
  struct constant_scoring: differentiable_scoring<queue_t> {
    dynamic_molecule protein;
    std::vector<int> active_updates;
    std::vector<int> initial_selection;
    bool capture_selection = true;
    bool coordinate_penalty = false;
    int atoms_per_ligand = 0;

    explicit constant_scoring(std::shared_ptr<scratchpad<queue_t>> scratch)
        : differentiable_scoring<queue_t>(scratch) {
      protein.resize(1, 0);
      protein.x(0) = protein.y(0) = protein.z(0) = fp_type{0};
    }
    dynamic_molecule& get_protein() { return protein; }
    void prepare(batch<static_molecule>& batch) override {
      atoms_per_ligand = batch.batch_max_atoms;
      const auto count = batch.num_ligands * this->scratch->configuration.population_number;
      initialize_buffer<buffer_data_type::SCORES>(this->scratch, count, fp_type{10}, true);
      active_updates.assign(count, 0);
    }
    void operator()() override {
      if (capture_selection) {
        auto& selected = this->scratch->template get<buffer_data_type::LS_SELECTED>();
        initial_selection.assign(selected(), selected() + selected.num_elements());
        capture_selection = false;
      }
      auto& scores = this->scratch->template get<buffer_data_type::SCORES>();
      std::fill_n(scores(), scores.num_elements(), fp_type{10});
      if (coordinate_penalty) {
        auto& x = this->scratch->template get<buffer_data_type::X_SCRATCH>();
        for (std::size_t i = 0; i < scores.num_elements(); ++i) {
          const fp_type displacement = x()[i * static_cast<std::size_t>(atoms_per_ligand)];
          scores()[i] += displacement * displacement;
        }
      }
    }
    void compute_gradient() override {
      auto& gradients = this->scratch->template get<buffer_data_type::GRADIENTS>();
      auto& selected = this->scratch->template get<buffer_data_type::LS_SELECTED>();
      for (std::size_t i = 0; i < gradients.num_elements(); ++i) {
        gradients()[i] = chromosome{};
        gradients()[i][0] = fp_type{1};
        active_updates[i] += selected()[i] != 0;
      }
    }
    void teardown_impl(batch<static_molecule>&) override {}
  };

  struct search_fixture {
    std::shared_ptr<scratchpad<queue_cpp>> scratch;
    std::shared_ptr<constant_scoring<queue_cpp>> score;
    batch<static_molecule> ligands;
    std::unique_ptr<adadelta<queue_cpp, constant_scoring>> search;

    explicit search_fixture(const knobs& conf)
        : scratch(std::make_shared<scratchpad<queue_cpp>>(conf, 0, device_type::CPU)),
          score(std::make_shared<constant_scoring<queue_cpp>>(scratch)) {
      ligands.num_ligands = 1;
      ligands.batch_max_atoms = 4;
      ligands.molecules[0] = std::make_unique<static_molecule>();
      auto& ligand = *ligands.molecules[0];
      ligand.resize(4, 3);
      ligand.properties.assign(property_type::NAME, "convergence_test");
      for (int atom = 0; atom < 4; ++atom) {
        ligand.x(atom) = static_cast<fp_type>(atom);
        ligand.y(atom) = ligand.z(atom) = fp_type{0};
        if (atom < 3) ligand.bonds(atom) = bond{atom, atom + 1};
      }
      initialize_buffer<buffer_data_type::CHROMOSOMES>(scratch, conf.population_number, chromosome{});
      search = std::make_unique<adadelta<queue_cpp, constant_scoring>>(
          scratch, score, ligand_placement{ligand_placement_mode::preserve, std::nullopt});
      search->prepare(ligands);
    }

    std::vector<int> selected() {
      auto& buffer = scratch->get<buffer_data_type::LS_SELECTED>();
      return {buffer(), buffer() + buffer.num_elements()};
    }
  };

  void repeated_searches_reactivate_individuals() {
    for (const std::size_t population : {std::size_t{1}, std::size_t{3}}) {
      knobs conf;
      conf.population_number = population;
      conf.num_generations = population == 1 ? 1 : 2;
      conf.lsit = MIN_ITERATIONS_TO_WARM_UP + MAX_BEST_SCORE_AGE + 2;
      conf.seed = 123;
      search_fixture f(conf);
      for (int invocation = 1; invocation <= 2; ++invocation) {
        (*f.search)();
        for (const auto updates : f.score->active_updates) {
          require(updates == invocation * MIN_ITERATIONS_TO_WARM_UP,
                  "Each local search call must reactivate stopped individuals and reset convergence");
        }
        const auto mask = f.selected();
        require(std::ranges::all_of(mask, [](const int active) { return active == 0; }),
                "Standalone and LGA searches must both apply individual autostop");
        auto& chromosomes = f.scratch->get<buffer_data_type::CHROMOSOMES>();
        require(std::ranges::all_of(std::span(chromosomes(), chromosomes.num_elements()),
                                   [](const chromosome& genes) { return genes == chromosome{}; }),
                "Search must return the best pose after a flat score trajectory");
      }
    }
  }

  void reproducible_generation_selection() {
    knobs conf;
    conf.population_number = 64;
    conf.num_generations = 2;
    conf.lsrate = 50;
    conf.lsit = 0;
    conf.seed = 123;
    search_fixture first(conf);
    search_fixture second(conf);
    (*first.search)();
    (*second.search)();
    const auto first_generation = first.score->initial_selection;
    require(first_generation == second.score->initial_selection,
            "Fixed seed must reproduce local search selection");
    first.score->capture_selection = second.score->capture_selection = true;
    (*first.search)();
    (*second.search)();
    require(first.score->initial_selection == second.score->initial_selection,
            "Fixed seed must reproduce selection across generations");
    require(first.score->initial_selection != first_generation,
            "Random selection must be resampled for each generation");
  }

  void budget_limit_refreshes_pose_and_score() {
    knobs conf;
    conf.population_number = 3;
    conf.num_generations = 2;
    conf.lsit = 2;
    search_fixture f(conf);
    f.score->coordinate_penalty = true;
    (*f.search)();
    auto& chromosomes = f.scratch->get<buffer_data_type::CHROMOSOMES>();
    auto& coordinates = f.scratch->get<buffer_data_type::X_SCRATCH>();
    auto& scores = f.scratch->get<buffer_data_type::SCORES>();
    for (std::size_t individual = 0; individual < conf.population_number; ++individual) {
      require(chromosomes()[individual] == chromosome{}, "Budget limit must restore the best chromosome");
      require(scores()[individual] == fp_type{10}, "Returned score must describe the restored pose");
      for (int atom = 0; atom < f.ligands.batch_max_atoms; ++atom) {
        const auto index = individual * static_cast<std::size_t>(f.ligands.batch_max_atoms) +
                           static_cast<std::size_t>(atom);
        require(coordinates()[index] == f.ligands.molecules[0]->x(atom),
                "Returned coordinates must describe the restored chromosome");
      }
    }
  }

  void best_individual_selection_is_refreshed() {
    knobs conf;
    conf.population_number = 4;
    conf.num_generations = 2;
    conf.ls_on_best = true;
    conf.lsrate = 50;
    conf.lsit = 0;
    search_fixture f(conf);
    auto& scores = f.scratch->get<buffer_data_type::SCORES>();
    const std::vector<fp_type> first_scores{4, 1, 3, 2};
    std::copy(first_scores.begin(), first_scores.end(), scores());
    (*f.search)();
    require(f.score->initial_selection == std::vector<int>({0, 1, 0, 1}),
            "Best-only local search must select the best individuals at the requested rate");
    const std::vector<fp_type> next_scores{1, 4, 2, 3};
    std::copy(next_scores.begin(), next_scores.end(), scores());
    f.score->capture_selection = true;
    (*f.search)();
    require(f.score->initial_selection == std::vector<int>({1, 0, 1, 0}),
            "Best-only selection must use the current generation's scores");
  }
} // namespace

int main() {
  try {
    warmup_and_best_pose();
    patience_and_improvement();
    individual_and_ligand_isolation();
    final_score_and_nonfinite_scores();
    stage_reset_and_pingpong();
    repeated_searches_reactivate_individuals();
    reproducible_generation_selection();
    budget_limit_refreshes_pose_and_score();
    best_individual_selection_is_refreshed();
  } catch (const std::exception& error) {
    std::cerr << "AdaDelta convergence regression failed: " << error.what() << '\n';
    return 1;
  }
  std::cout << "AdaDelta convergence regressions passed\n";
}
