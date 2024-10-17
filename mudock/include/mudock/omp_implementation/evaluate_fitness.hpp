#pragma once

#include <cstddef>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/grid.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  void evaluate_fitness(const int batch_ligands,
                        const int num_generations,
                        const int tournament_length,
                        const fp_type mutation_prob,
                        const int chromosome_number,
                        const int chromosome_stride,
                        const int atom_stride,
                        const int rotamers_stride,
                        const int nonbond_stride,
                        const fp_type* __restrict__ original_ligand_x,
                        const fp_type* __restrict__ original_ligand_y,
                        const fp_type* __restrict__ original_ligand_z,
                        fp_type* __restrict__ scratch_ligand_x,
                        fp_type* __restrict__ scratch_ligand_y,
                        fp_type* __restrict__ scratch_ligand_z,
                        const fp_type* __restrict__ ligand_vol,
                        const fp_type* __restrict__ ligand_solpar,
                        const fp_type* __restrict__ ligand_charge,
                        const int* __restrict__ ligand_num_hbond,
                        const fp_type* __restrict__ ligand_Rij_hb,
                        const fp_type* __restrict__ ligand_Rii,
                        const fp_type* __restrict__ ligand_epsij_hb,
                        const fp_type* __restrict__ ligand_epsii,
                        const int* __restrict__ ligand_num_nonbonds,
                        const int* __restrict__ ligand_nonbond_a1,
                        const int* __restrict__ ligand_nonbond_a2,
                        const int* __restrict__ ligand_num_atoms,
                        const int* __restrict__ ligand_num_rotamers,
                        const int* __restrict__ ligand_fragments,
                        const int* __restrict__ frag_start_atom_index,
                        const int* __restrict__ frag_stop_atom_index,
                        chromosome* __restrict__ chromosomes,
                        // TODO input data maps
                        // const cudaTextureObject_t* atom_textures,
                        // const int* atom_tex_indexes,
                        // const cudaTextureObject_t electro_texture,
                        // const cudaTextureObject_t desolv_texture,
                        // TODO random state
                        fp_type* ligand_scores,
                        chromosome* best_chromosomes);
} // namespace mudock
