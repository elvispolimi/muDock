#include <cstdio>
#include <mudock/cuda_implementation/mutate.cuh>
#include <mudock/cuda_implementation/evaluate_fitness.cuh>

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
                                   fp_type* ligand_scores) {
    const int ligand_id = blockIdx.x;

    const std::size_t num_atoms    = ligand_num_atoms[ligand_id];
    const std::size_t num_rotamers = ligand_num_rotamers[ligand_id];

    fp_type* l_ligand_x           = ligand_x + population_stride * ligand_id * atom_stride;
    fp_type* l_ligand_y           = ligand_y + population_stride * ligand_id * atom_stride;
    fp_type* l_ligand_z           = ligand_z + population_stride * ligand_id * atom_stride;
    const individual* l_population      = population + population_stride * ligand_id;
    const auto* l_fragments             = ligand_fragments + ligand_id * atom_stride;
    const auto* l_frag_start_atom_index = frag_start_atom_index + ligand_id * rotamers_stride;
    const auto* l_frag_stop_atom_index  = frag_stop_atom_index + ligand_id * rotamers_stride;

    for (int el = 0; el < population_stride; ++el) {
      const individual* l_individual = l_population+el;
      // Modify coordinates
      apply_cuda(l_ligand_x,
            l_ligand_y,
            l_ligand_z,
            l_individual->genes.data(),
            l_fragments,
            l_frag_start_atom_index,
            l_frag_stop_atom_index,
            num_rotamers,
            atom_stride,
            num_atoms);
    }
  }
} // namespace mudock
