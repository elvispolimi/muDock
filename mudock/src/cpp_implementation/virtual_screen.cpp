#include <mudock/cpp_implementation/virtual_screen.hpp>

namespace mudock {

  virtual_screen_cpp::virtual_screen_cpp([[maybe_unused]] std::shared_ptr<dynamic_molecule>& protein)
      : generator(protein->num_atoms()), dist(fp_type{0}, fp_type{10}) {
    // NOTE: at the moment we don't care about the protein
  }

  void virtual_screen_cpp::operator()(static_molecule& ligand) {
    ligand.properties.assign(property_type::SCORE, std::to_string(dist(generator)));
  }

} // namespace mudock
