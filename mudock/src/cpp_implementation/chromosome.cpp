#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>

namespace mudock {

  void apply(std::span<fp_type> x,
             std::span<fp_type> y,
             std::span<fp_type> z,
             const chromosome& c,
             const fragments<static_containers>& fragments) {
    // apply rigid transformations
    translate_molecule(x, y, z, c[0], c[1], c[2]);
    rotate_molecule(x, y, z, c[3], c[4], c[5]);

    // change the molecule shape
    const auto num_rotamers = fragments.get_num_rotatable_bonds();
    for (std::size_t i = 0; i < num_rotamers; ++i) {
      const auto bitmask                   = fragments.get_mask(i);
      const auto [start_index, stop_index] = fragments.get_rotatable_atoms(i);
      rotate_fragment(x, y, z, bitmask, start_index, stop_index, c[std::size_t{6} + i]);
    }
  }
} // namespace mudock
