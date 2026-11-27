#pragma once

#include <mudock/type_alias.hpp>
#include <type_traits>

namespace mudock {
  struct adt_score_args {
    int batch_ligands;
    int scores_per_ligand;
    int ligand_index;
    int batch_atom_stride;
    int map_index_x;
    int map_index_xy;
    int map_index_xyz;
    const fp_type* vols;
    const fp_type* solpars;
    const fp_type* charges;
    const int* num_atoms;
    const int* num_rotamers;
    const int* num_nonbonds;
    const int* nonbond_a1;
    const int* nonbond_a2;
    const fp_type* nonbond_cA;
    const fp_type* nonbond_cB;
    const int* nonbond_xB;
    const fp_type* grid_maps;
    const int* atom_map_offsets;
    const fp_type* minimum;
    const fp_type* maximum;
    const fp_type* center;
    fp_type* scores;
    const fp_type* x;
    const fp_type* y;
    const fp_type* z;
  };
  static_assert(std::is_trivially_copyable_v<adt_score_args>);
} // namespace mudock
