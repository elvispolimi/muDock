#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/g_highway_implementation/geometric_transformations.hpp>
#include <mudock/molecule/fragments.hpp>

namespace mudock {
  void apply(fp_type* __restrict__ x,
             fp_type* __restrict__ y,
             fp_type* __restrict__ z,
             const chromosome& c,
             const int num_atoms,
             const int num_rotamers,
             const int* __restrict__ frag_masks,
             const int* __restrict__ frag_start_indexes,
             const int* __restrict__ frag_stop_indexes) {
    // apply rigid transformations
    translate_molecule(x, y, z, num_atoms, c[0], c[1], c[2]);
    rotate_molecule(x, y, z, num_atoms, c[3], c[4], c[5]);

    for (int i = 0; i < num_rotamers; ++i) {
      const auto* bitmask    = frag_masks + i * num_atoms;
      const auto start_index = frag_start_indexes[i];
      const auto stop_index  = frag_stop_indexes[i];
      rotate_fragment(x, y, z, num_atoms, bitmask, start_index, stop_index, c[int{6} + i]);
    }
  }

} // namespace mudock
