#include <mudock/cpp_kernels/virtual_screen.hpp>

namespace mudock {

  virtual_screen_cpp::virtual_screen_cpp(std::shared_ptr<dynamic_molecule>& protein)
      : kernel_interface(protein), generator(protein->num_atoms()), dist(fp_type{0}, fp_type{10}) {}

  void virtual_screen_cpp::virtual_screen(const std::span<std::unique_ptr<static_molecule>> ligands) {
    for (auto& ligand_ptr: ligands) {
      ligand_ptr->properties.assign(property_type::SCORE, std::to_string(dist(generator)));
    }
  }

} // namespace mudock
