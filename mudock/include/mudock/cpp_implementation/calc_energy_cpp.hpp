#pragma once

#include <mudock/chem/ligand_maps.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/type_alias.hpp>
#include <stdint.h>

namespace mudock {
  // TODO template parameters on num_atoms?
  // Buckets?
  void evaluate_fitness(const fp_type* __restrict__ ligand_x,
                        const fp_type* __restrict__ ligand_y,
                        const fp_type* __restrict__ ligand_z,
                        const fp_type* __restrict__ ligand_vol,
                        const fp_type* __restrict__ ligand_solpar,
                        const fp_type* __restrict__ ligand_charge,
                        const int* __restrict__ ligand_num_hbond,
                        const fp_type* __restrict__ ligand_Rij_hb,
                        const fp_type* __restrict__ ligand_Rii,
                        const fp_type* __restrict__ ligand_epsij_hb,
                        const fp_type* __restrict__ ligand_epsii,
                        const ligand_map_types* __restrict__ map_ligand_types,
                        const int num_atoms,
                        const int num_rotamers,
                        const int* __restrict__ frag_masks,
                        const int* __restrict__ frag_start_indexes,
                        const int* __restrict__ frag_stop_indexes,
                        const uint_fast8_t* __restrict__ nbmatrix,
                        const fp_type* const __restrict__* const __restrict__ grid_maps,
                        const fp_type* __restrict__ electro_map,
                        const fp_type* __restrict__ desolv_map,
                        const int num_generations,
                        const int population_size,
                        const int tournament_length,
                        const fp_type mutation_prob,
                        const fp_type* __restrict__ minimum,
                        const fp_type* __restrict__ maximum,
                        const fp_type* __restrict__ center,
                        const int map_index_x,
                        const int map_index_xy,
                        individual* __restrict__ population_buffer1,
                        individual* __restrict__ population_buffer2,
                        const int seed);
} // namespace mudock
