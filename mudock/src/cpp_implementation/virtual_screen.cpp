#include <algorithm>
#include <mudock/cpp_implementation/calc_energy_cpp.hpp>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/trilinear_interpolation.hpp>
#include <mudock/cpp_implementation/virtual_screen.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/grid.hpp>
#include <mudock/molecule.hpp>
#include <mudock/utils.hpp>
#include <vector>

namespace mudock {

  virtual_screen_cpp::virtual_screen_cpp(std::shared_ptr<const grid_atom_mapper> &_grid_atom_maps,
                                         std::shared_ptr<const grid_map> &_electro_map,
                                         std::shared_ptr<const grid_map> &_desolv_map,
                                         const knobs &knobs)
      : grid_atom_maps(_grid_atom_maps),
        electro_map(_electro_map),
        desolv_map(_desolv_map),
        population(knobs.population_number),
        next_population(knobs.population_number),
        configuration(knobs) {}

  template<typename T>
  [[nodiscard]] const T virtual_screen_cpp::random_gen_cpp(const T &min, const T &max) {
    fp_type value;
    if constexpr (is_debug())
      // TODO value here for debug
      value = fp_type{0.4};
    else {
      value = dist(generator);
    }
    return static_cast<T>(value * (max - min) + min);
  }

  int virtual_screen_cpp::get_selection_distribution() {
    return random_gen_cpp<int>(0, configuration.population_number - 1);
  };

  fp_type virtual_screen_cpp::get_init_change_distribution() { return random_gen_cpp<fp_type>(-45, 45); }
  fp_type virtual_screen_cpp::get_mutation_change_distribution() { return random_gen_cpp<fp_type>(-10, 10); };
  fp_type virtual_screen_cpp::get_mutation_coin_distribution() { return random_gen_cpp<fp_type>(0, 1); };
  int virtual_screen_cpp::get_crossover_distribution(const int &num_rotamers) {
    return random_gen_cpp<int>(0, 6 + num_rotamers);
  };

  const chromosome &virtual_screen_cpp::tournament_selection() {
    const auto num_iterations = configuration.tournament_length;
    auto best_individual      = get_selection_distribution();
    for (std::size_t i = 0; i < num_iterations; ++i) {
      auto contendent = get_selection_distribution();
      if (population[contendent].score < population[best_individual].score) {
        best_individual = contendent;
      }
    }
    return population[best_individual].genes;
  }

  void virtual_screen_cpp::operator()(static_molecule &ligand) {
    // Reset the random number generator to improve consistency
    generator = std::mt19937{static_cast<size_t>(ligand.num_atoms())};

    // Place the molecule to the center of the target protein
    const auto x = ligand.get_x(), y = ligand.get_y(), z = ligand.get_z();
    const auto ligand_center_of_mass = compute_center_of_mass(x, y, z);
    translate_molecule(x,
                       y,
                       z,
                       electro_map->center.x - ligand_center_of_mass.x,
                       electro_map->center.y - ligand_center_of_mass.y,
                       electro_map->center.z - ligand_center_of_mass.z);

    // Find out the rotatable bonds in the ligand
    auto graph = make_graph(ligand.get_bonds());
    const fragments<static_containers> ligand_fragments{graph, ligand.get_bonds(), ligand.num_atoms()};

    const auto coordinate_step = fp_type{0.2};
    const auto angle_step      = fp_type{4};

    // Randomly initialize the population
    const auto num_rotamers = ligand_fragments.get_num_rotatable_bonds();
    for (auto &element: population) {
      for (int i{0}; i < 3; ++i) { // initialize the rigid translation
        element.genes[i] = get_init_change_distribution() * coordinate_step;
      }
      for (int i{3}; i < 6 + num_rotamers; ++i) { // initialize the rotations
        element.genes[i] = get_init_change_distribution() * angle_step;
      }
    }

    // Get weed bonds and non bonds lists
    const int num_atoms = ligand.num_atoms();
    grid<uint_fast8_t, index2D> nbmatrix{{num_atoms, num_atoms}};
    nonbonds(nbmatrix, ligand.get_bonds(), num_atoms);
    std::vector<non_bond_parameter> non_bond_list;
    weed_bonds(nbmatrix, non_bond_list, num_atoms, ligand_fragments);

    // Simulate the population evolution for the given amount of time
    const auto num_generations = configuration.num_generations;
    for (std::size_t generation = 0; generation < num_generations; ++generation) {
      // Evaluate the fitness of the population
      for (std::size_t element_index = 0; element_index < population.size(); ++element_index) {
        auto &element = population[element_index];
        // copy the ligand original coordinates in temporary array
        auto altered_x = static_containers::atoms_size<fp_type>{};
        auto altered_y = static_containers::atoms_size<fp_type>{};
        auto altered_z = static_containers::atoms_size<fp_type>{};
        std::copy(std::cbegin(x), std::cend(x), std::begin(altered_x));
        std::copy(std::cbegin(y), std::cend(y), std::begin(altered_y));
        std::copy(std::cbegin(z), std::cend(z), std::begin(altered_z));
        // TODO check it it makes sense -> print the MOL2
        // apply the transformation encoded in the element genes to the original ligand
        apply(std::span(std::begin(altered_x), x.size()),
              std::span(std::begin(altered_y), y.size()),
              std::span(std::begin(altered_z), z.size()),
              element.genes,
              ligand_fragments);

        // compute the energy of the system
        // printf("%ld %ld ", generation, element_index);
        const auto energy = calc_energy(altered_x,
                                        altered_y,
                                        altered_z,
                                        ligand.get_vol(),
                                        ligand.get_solpar(),
                                        ligand.get_charge(),
                                        ligand.get_num_hbond(),
                                        ligand.get_Rij_hb(),
                                        ligand.get_Rii(),
                                        ligand.get_epsij_hb(),
                                        ligand.get_epsii(),
                                        ligand.get_autodock_type(),
                                        ligand.num_atoms(),
                                        ligand_fragments.get_num_rotatable_bonds(),
                                        non_bond_list,
                                        *grid_atom_maps,
                                        *electro_map,
                                        *desolv_map);
        element.score     = energy; // dummy implementation to test the genetic
      }

      // Generate the new population
      for (auto &next_individual: next_population) {
        // select the parent
        const auto &parent1 = tournament_selection();
        const auto &parent2 = tournament_selection();

        // generate the offspring
        const auto split_index = get_crossover_distribution(num_rotamers);
        std::copy(std::begin(parent1), std::begin(parent1) + split_index, std::begin(next_individual.genes));
        std::copy(std::begin(parent2) + split_index,
                  std::end(parent2),
                  std::begin(next_individual.genes) + split_index);
        next_individual.score = fp_type{0};

        // mutate the offspring
        for (int i{0}; i < 3; ++i) {
          if (get_mutation_coin_distribution() < configuration.mutation_prob)
            next_individual.genes[i] += get_mutation_change_distribution() * coordinate_step;
        }
        for (int i{3}; i < 6 + num_rotamers; ++i) {
          if (get_mutation_coin_distribution() < configuration.mutation_prob)
            next_individual.genes[i] += get_mutation_change_distribution() * angle_step;
        }
      }

      // swap the new population with the old one
      population.swap(next_population);
    }

    // update the ligand position with the best one that we found
    const auto best_individual_it =
        std::min_element(std::begin(next_population),
                         std::end(next_population),
                         [](const auto a, const auto b) { return a.score < b.score; });
    apply(x, y, z, best_individual_it->genes, ligand_fragments);
    ligand.properties.assign(property_type::SCORE, std::to_string(best_individual_it->score));
  }

} // namespace mudock
