#pragma once

#include <array>
#include <cuda.h>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/molecule.hpp>
#include <mudock/omp_implementation/geometric_transformations.hpp>
#include <mudock/type_alias.hpp>
#include <random>
#include <span>

namespace mudock {
  void inline apply_omp(fp_type* __restrict__ x,
                        fp_type* __restrict__ y,
                        fp_type* __restrict__ z,
                        const chromosome& chromosome,
                        const int* __restrict__ fragments,
                        const int* __restrict__ fragments_start_index,
                        const int* __restrict__ fragments_stop_index,
                        const int num_rotamers,
                        const int stride_atoms,
                        const int num_atoms) {
    fp_type reference_x{0}, reference_y{0}, reference_z{0};
#pragma omp parallel for reduction(+ : reference_x, reference_y, reference_z)
    for (int atom_index = 0; atom_index < num_atoms; ++atom_index) {
      reference_x += x[atom_index];
      reference_y += y[atom_index];
      reference_z += z[atom_index];
    }
    reference_x /= num_atoms;
    reference_y /= num_atoms;
    reference_z /= num_atoms;

    // change the molecule shape
    for (int i = 0; i < num_rotamers; ++i) {
      const int* bitmask = fragments + i * stride_atoms;
      rotate_fragment_omp(x,
                          y,
                          z,
                          bitmask,
                          fragments_start_index[i],
                          fragments_stop_index[i],
                          &chromosome[6 + i],
                          num_atoms);
    }

    recenter_molecule_omp(x, y, z, num_atoms, reference_x, reference_y, reference_z);

    // Apply rigid transformations after the conformation is recentered.
    rotate_molecule_omp(x, y, z, &chromosome[3], &chromosome[4], &chromosome[5], num_atoms);
    translate_molecule_omp(x, y, z, &chromosome[0], &chromosome[1], &chromosome[2], num_atoms);
  }

} // namespace mudock
