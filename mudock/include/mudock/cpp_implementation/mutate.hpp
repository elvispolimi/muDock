#pragma once

#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<cpu_vectorization vect>
  void translate_molecule(fp_type* __restrict__ x,
                          fp_type* __restrict__ y,
                          fp_type* __restrict__ z,
                          const int num_atoms,
                          const fp_type offset_x,
                          const fp_type offset_y,
                          const fp_type offset_z);

  template<cpu_vectorization vect>
  void rotate_molecule(fp_type* __restrict__ x,
                       fp_type* __restrict__ y,
                       fp_type* __restrict__ z,
                       const int num_atoms,
                       const fp_type angle_x,
                       const fp_type angle_y,
                       const fp_type angle_z);

  template<cpu_vectorization vect>
  void rotate_fragment(fp_type* __restrict__ x,
                       fp_type* __restrict__ y,
                       fp_type* __restrict__ z,
                       const int num_atoms,
                       const int* __restrict__ frag_mask,
                       const int start_index,
                       const int stop_index,
                       const fp_type angle);

  template<cpu_vectorization vect>
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
    translate_molecule<vect>(x, y, z, num_atoms, c[0], c[1], c[2]);
    rotate_molecule<vect>(x, y, z, num_atoms, c[3], c[4], c[5]);

    // change the molecule shape
    for (int i = 0; i < num_rotamers; ++i) {
      const auto* bitmask    = frag_masks + i * num_atoms;
      const auto start_index = frag_start_indexes[i];
      const auto stop_index  = frag_stop_indexes[i];
      rotate_fragment<vect>(x, y, z, num_atoms, bitmask, start_index, stop_index, c[int{6} + i]);
    }
  }
} // namespace mudock
