#include <mudock/cpp_implementation/virtual_screen.hpp>
#include <mudock/kernels.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  virtual_screen_cpp::virtual_screen_cpp([[maybe_unused]] std::shared_ptr<const dynamic_molecule> &_protein,
                                         std::shared_ptr<const grid_atom_mapper> _grid_atom_maps,
                                         std::shared_ptr<const grid_map> _electro_map,
                                         std::shared_ptr<const grid_map> _desolv_map,
                                         const knobs knobs)
      : protein(_protein),
        grid_atom_maps(_grid_atom_maps),
        electro_map(_electro_map),
        desolv_map(_desolv_map),
        generator(_protein->num_atoms()),
        configuration(knobs) {
    // NOTE: at the moment we don't care about the protein
  }

  void virtual_screen_cpp::operator()(static_molecule &ligand) {
    // TODO tentative fragments creation
    // Should it be like this?
    // Graph creation is done two times actually
    // Assign fragments once in convert function?
    auto graph = make_graph(ligand.get_bonds());
    const fragments<static_containers> ligand_fragments{graph, ligand.get_bonds(), ligand.num_atoms()};

    const auto energy =
        calc_energy(*protein, ligand, ligand_fragments, *grid_atom_maps, *electro_map, *desolv_map);

    ligand.properties.assign(property_type::SCORE, std::to_string(energy));
  }

} // namespace mudock
