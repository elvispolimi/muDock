#pragma once

#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/type_alias.hpp>
#include <type_traits>

namespace mudock {
  struct geom_transform_args {
    int batch_ligands;
    int batch_atom_stride;
    int chromosomes_per_ligand;
    const chromosome* chromosomes;
    const fp_type* original_x;
    const fp_type* original_y;
    const fp_type* original_z;
    const int* num_atoms;
    const int* num_rotamers;
    const int* fragments;
    const int* ligand_fragments_start;
    const int* fragments_start_index;
    const int* fragments_stop_index;
    const int* frag_indices_start;
    fp_type* x_scratch;
    fp_type* y_scratch;
    fp_type* z_scratch;
  };
  static_assert(std::is_trivially_copyable_v<geom_transform_args>);
} // namespace mudock
