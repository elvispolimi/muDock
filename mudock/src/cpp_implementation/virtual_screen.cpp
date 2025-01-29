#include <memory>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/cpp_implementation/calc_energy_cpp.hpp>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/virtual_screen.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/grid.hpp>
#include <mudock/molecule.hpp>
#include <mudock/scorep_utils.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  virtual_screen_cpp::virtual_screen_cpp(std::shared_ptr<const grid_atom_mapper>& _grid_atom_maps,
                                         std::shared_ptr<const grid_map>& _electro_map,
                                         std::shared_ptr<const grid_map>& _desolv_map,
                                         const knobs& knobs)
      : grid_atom_maps(_grid_atom_maps),
        electro_map(_electro_map),
        desolv_map(_desolv_map),
        population(knobs.population_number),
        next_population(knobs.population_number),
        configuration(knobs) {}

  void virtual_screen_cpp::operator()(static_molecule& ligand) {
    SCOREP_MARKER_START(ga, "GA");
    const auto seed =
        static_cast<size_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    // Place the molecule to the center of the target protein
    const int num_atoms = ligand.num_atoms();
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
    const auto ligand_fragments =
        std::make_unique<fragments<static_containers>>(graph, ligand.get_bonds(), ligand.num_atoms());

    const auto num_rotamers = ligand_fragments.get()->get_num_rotatable_bonds();

    // Get weed bonds and non bonds lists
    grid<uint_fast8_t, index2D> nbmatrix{{num_atoms, num_atoms}};
    nonbonds(nbmatrix, ligand.get_bonds(), num_atoms);
    std::vector<int> non_bond_list_a1, non_bond_list_a2;
    weed_bonds(nbmatrix, non_bond_list_a1, non_bond_list_a2, num_atoms, *ligand_fragments.get());
    // weed_bonds(nbmatrix, num_atoms, *ligand_fragments.get());

    const auto non_bond_size = non_bond_list_a1.size();
    std::vector<fp_type> cA_v, cB_v;
    std::vector<int> xB_v;
    cA_v.resize(non_bond_size);
    cB_v.resize(non_bond_size);
    xB_v.resize(non_bond_size);
    for(int index=0; index<non_bond_size; ++index){
      const int& a1 = non_bond_list_a1[index];
      const int& a2 = non_bond_list_a2[index];

      const auto& hbond_i    = ligand.num_hbond(a1);
      const auto& hbond_j    = ligand.num_hbond(a2);
      const auto& Rij_hb_i   = ligand.Rij_hb(a1);
      const auto& Rij_hb_j   = ligand.Rij_hb(a2);
      const auto& Rii_i      = ligand.Rii(a1);
      const auto& Rii_j      = ligand.Rii(a2);
      const auto& epsij_hb_i = ligand.epsij_hb(a1);
      const auto& epsij_hb_j = ligand.epsij_hb(a2);
      const auto& epsii_i    = ligand.epsii(a1);
      const auto& epsii_j    = ligand.epsii(a2);

      // we need to determine the correct xA and xB exponents
      const int xA = xA_default; // for both LJ, 12-6 and HB, 12-10, xA is 12
      int xB       = xB_default;  // assume we have LJ, 12-6

      fp_type Rij{(Rii_i + Rii_j) * fp_type{0.5}}, epsij{std::sqrt(epsii_i * epsii_j)};
      if ((hbond_i == 1 || hbond_i == 2) && hbond_j > 2) {
        // i is a donor and j is an acceptor.
        // i is a hydrogen, j is a heteroatom
        Rij   = Rij_hb_j;
        epsij = epsij_hb_j;
        xB    = 10;
      } else if ((hbond_i > 2) && (hbond_j == 1 || hbond_j == 2)) {
        // i is an acceptor and j is a donor.
        // i is a heteroatom, j is a hydrogen
        Rij   = Rij_hb_i;
        epsij = epsij_hb_i;
        xB    = 10;
      }
      fp_type cA{0};
      fp_type cB{0};
      if (xA != xB) {
        const fp_type tmp = epsij / (xA - xB);
        cA = tmp * std::pow(Rij, static_cast<fp_type>(xA)) * xB;
        cB = tmp * std::pow(Rij, static_cast<fp_type>(xB)) * xA;
      }
      cA_v[index] = cA;
      cB_v[index] = cB;
      xB_v[index] = xB;
    }

    const fp_type minimum[3] = {electro_map.get()->minimum_coord.x,
                                electro_map.get()->minimum_coord.y,
                                electro_map.get()->minimum_coord.z};
    const fp_type maximum[3] = {electro_map.get()->maximum_coord.x,
                                electro_map.get()->maximum_coord.y,
                                electro_map.get()->maximum_coord.z};
    const fp_type center[3]  = {electro_map.get()->center.x,
                                electro_map.get()->center.y,
                                electro_map.get()->center.z};

    fp_type const* grid_maps[num_ligand_map_types()];
    constexpr_for<0, num_ligand_map_types(), 1>([&](const int type) {
      grid_maps[type] = grid_atom_maps.get()
                            ->get_atom_map(autodock_type_from_map(static_cast<ligand_map_types>(type)))
                            .data();
    });

    std::vector<int> frag_masks;
    std::vector<int> frag_start_indexes;
    std::vector<int> frag_stop_indexes;
    frag_masks.resize(num_atoms * num_rotamers);
    frag_start_indexes.resize(num_rotamers);
    frag_stop_indexes.resize(num_rotamers);
    for (int rot = 0; rot < num_rotamers; ++rot) {
      std::memcpy((frag_masks.data() + num_atoms * rot),
                  ligand_fragments.get()->get_mask(rot).data(),
                  num_atoms * sizeof(int));
      const auto [start_index, stop_index] = ligand_fragments.get()->get_rotatable_atoms(rot);
      frag_start_indexes.data()[rot]       = start_index;
      frag_stop_indexes.data()[rot]        = stop_index;
    }

    std::vector<int> map_ligand_types;
    map_ligand_types.resize(num_atoms);
    for (int i = 0; i < num_atoms; i++)
      map_ligand_types[i] = static_cast<int>(map_from_autodock_type(ligand.autodock_type(i)));

    // Simulate the population evolution for the given amount of time
    evaluate_fitness(x.data(),
                     y.data(),
                     z.data(),
                     ligand.get_vol().data(),
                     ligand.get_solpar().data(),
                     ligand.get_charge().data(),
                    //  ligand.get_num_hbond().data(),
                    //  ligand.get_Rij_hb().data(),
                    //  ligand.get_Rii().data(),
                    //  ligand.get_epsij_hb().data(),
                    //  ligand.get_epsii().data(),
                     map_ligand_types.data(),
                     num_atoms,
                     ligand_fragments.get()->get_num_rotatable_bonds(),
                     frag_masks.data(),
                     frag_start_indexes.data(),
                     frag_stop_indexes.data(),
                     non_bond_size,
                     non_bond_list_a1.data(),
                     non_bond_list_a2.data(),
                     cA_v.data(),
                     cB_v.data(),
                     xB_v.data(),
                     grid_maps,
                     electro_map.get()->data(),
                     desolv_map.get()->data(),
                     configuration.num_generations,
                     configuration.population_number,
                     configuration.tournament_length,
                     configuration.mutation_prob,
                     minimum,
                     maximum,
                     center,
                     electro_map.get()->index.size_x(),
                     electro_map.get()->index.size_xy(),
                     population.data(),
                     next_population.data(),
                     seed);

    // update the ligand position with the best one that we found
    const std::vector<individual>& last_population =
        (configuration.num_generations % 2 == 0) ? next_population : population;

    const auto best_individual_it =
        std::min_element(std::begin(last_population),
                         std::end(last_population),
                         [](const auto a, const auto b) { return a.score < b.score; });
    apply(x, y, z, best_individual_it->genes, *ligand_fragments.get());
    ligand.properties.assign(property_type::SCORE, std::to_string(best_individual_it->score));
    SCOREP_MARKER_STOP(ga);
  }

} // namespace mudock
