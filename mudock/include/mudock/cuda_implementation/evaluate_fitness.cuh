#pragma once

#include <cstddef>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/cuda_implementation/mutate.cuh>

namespace mudock {

  __global__ void evaluate_fitness(const std::size_t population_stride,
                                   const std::size_t atom_stride,
                                   const std::size_t rotamers_stride,
                                   fp_type* __restrict__ ligand_x,
                                   fp_type* __restrict__ ligand_y,
                                   fp_type* __restrict__ ligand_z,
                                   const fp_type* __restrict__ ligand_vol,
                                   const fp_type* __restrict__ ligand_solpar,
                                   const fp_type* __restrict__ ligand_charge,
                                   const std::size_t* __restrict__ ligand_num_hbond,
                                   const fp_type* __restrict__ ligand_Rij_hb,
                                   const fp_type* __restrict__ ligand_Rii,
                                   const fp_type* __restrict__ ligand_epsij_hb,
                                   const fp_type* __restrict__ ligand_epsii,
                                   const std::size_t* __restrict__ ligand_autodock_type,
                                   // TODO
                                   const fp_type* __restrict__ ligand_bond,
                                   const std::size_t* __restrict__ ligand_num_atoms,
                                   const std::size_t* __restrict__ ligand_num_rotamers,
                                   const int* __restrict__ ligand_fragments,
                                   const std::size_t* __restrict__ frag_start_atom_index,
                                   const std::size_t* __restrict__ frag_stop_atom_index,
                                   const individual* __restrict__ population,
                                   // TODO
                                   const fp_type* __restrict__ grid_maps,
                                   const fp_type* __restrict__ electro_map,
                                   const fp_type* __restrict__ desolv_map,
                                   fp_type* ligand_scores);
}
