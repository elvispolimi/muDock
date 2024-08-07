#include <mudock/kernels/calc_energy_cpp.hpp>

namespace mudock {

  void calc_energy(const static_molecule& receptor,
                   const dynamic_molecule& ligand,
                   const std::vector<grid_map> grid_maps) {
    for (size_t index = 0; index < ligand.num_atoms(); ++index) {
      const fp_type x = ligand.
    }
  }

} // namespace mudock
