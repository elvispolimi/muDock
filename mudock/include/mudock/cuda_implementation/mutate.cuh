#pragma once

#include <array>
#include <cuda.h>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cuda_implementation/geometric_transformations.cuh>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <random>
#include <span>

namespace mudock {
  template<int MAX_ATOMS, int MAX_ROTAMERS>
  __device__ void apply_cuda(fp_type* __restrict__ x,
                             fp_type* __restrict__ y,
                             fp_type* __restrict__ z,
                             const chromosome& chromosome,
                             const int* __restrict__ fragments,
                             const int* __restrict__ fragments_start_index,
                             const int* __restrict__ fragments_stop_index,
                             const int num_rotamers,
                             const int num_atoms) {
    // apply rigid transformations
    translate_molecule_cuda<MAX_ATOMS>(x, y, z, &chromosome[0], &chromosome[1], &chromosome[2], num_atoms);
    rotate_molecule_cuda<MAX_ATOMS>(x, y, z, &chromosome[3], &chromosome[4], &chromosome[5], num_atoms);
// #pragma unroll
//     for (int atom_index = threadIdx.x; atom_index < MAX_ATOMS; atom_index += blockDim.x) {
//       // for (int atom_index = local_thread_id; atom_index < num_atoms; atom_index += thread_per_block) {
//       if (atom_index < num_atoms) {
//         printf("CUDA Before %d %f %f %f\n", atom_index, x[atom_index], y[atom_index], z[atom_index]);
//       }
//     }

// change the molecule shape
#pragma unroll
    for (int i = 0; i < MAX_ROTAMERS; ++i) {
      if (i < num_rotamers) {
        const int* bitmask = fragments + i * num_atoms;
        // if (threadIdx.x == 0)
        //   printf("CUDA %d %d %d %d | %d %d\n",
        //          i,
        //          bitmask[0],
        //          bitmask[1],
        //          bitmask[2],
        //          fragments_start_index[i],
        //          fragments_stop_index[i]);
        rotate_fragment_cuda<MAX_ATOMS>(x,
                                        y,
                                        z,
                                        bitmask,
                                        fragments_start_index[i],
                                        fragments_stop_index[i],
                                        &chromosome[6 + i],
                                        num_atoms);
      }
    }
    // #pragma unroll
    //     for (int atom_index = threadIdx.x; atom_index < MAX_ATOMS; atom_index += blockDim.x) {
    //       // for (int atom_index = local_thread_id; atom_index < num_atoms; atom_index += thread_per_block) {
    //       if (atom_index < num_atoms) {
    //         printf("CUDA AFTER %d %f %f %f\n", atom_index, x[atom_index], y[atom_index], z[atom_index]);
    //       }
    //     }
  }

} // namespace mudock
