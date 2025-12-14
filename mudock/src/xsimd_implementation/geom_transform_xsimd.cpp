#include <cassert>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/xsimd_implementation/geom_transform_xsimd.hpp>
#include <mudock/xsimd_implementation/mutate_xsimd.hpp>

namespace mudock {
  template<>
  void geom_kernel<queue_xsimd>::operator()() {
    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      const int stride_atoms                     = ligand_index * batch_atoms;
      const int num_atoms                        = num_atoms_b[ligand_index];
      const int num_rotamers                     = num_rotamers_b[ligand_index];
      const fp_type* __restrict__ x_coords       = x_coords_b + stride_atoms;
      const fp_type* __restrict__ y_coords       = y_coords_b + stride_atoms;
      const fp_type* __restrict__ z_coords       = z_coords_b + stride_atoms;
      const chromosome* __restrict__ chromosomes = chromosomes_b + ligand_index * chromsomes_per_ligand;

      const int* fragments          = ligand_fragments_b + ligand_fragments_start_b[ligand_index];
      const int* frag_start_indices = frag_start_indices_b + frag_indices_start_b[ligand_index];
      const int* frag_stop_indices  = frag_stop_indices_b + frag_indices_start_b[ligand_index];

      fp_type* __restrict__ x_scratch = x_scratch_b + ligand_index * batch_atoms * chromsomes_per_ligand;
      fp_type* __restrict__ y_scratch = y_scratch_b + ligand_index * batch_atoms * chromsomes_per_ligand;
      fp_type* __restrict__ z_scratch = z_scratch_b + ligand_index * batch_atoms * chromsomes_per_ligand;
      for (int element_index = 0; element_index < chromsomes_per_ligand; ++element_index) {
        fp_type* __restrict__ x_scratch_chromosome = x_scratch + element_index * batch_atoms;
        fp_type* __restrict__ y_scratch_chromosome = y_scratch + element_index * batch_atoms;
        fp_type* __restrict__ z_scratch_chromosome = z_scratch + element_index * batch_atoms;

        std::copy(x_coords, x_coords + num_atoms, x_scratch_chromosome);
        std::copy(y_coords, y_coords + num_atoms, y_scratch_chromosome);
        std::copy(z_coords, z_coords + num_atoms, z_scratch_chromosome);

        // TODO try with loop perforation
        const chromosome& chrom = chromosomes[element_index];

        // FIXME with different vectorization type
        apply<cpu_vectorization::XSIMD>(x_scratch_chromosome,
                                        y_scratch_chromosome,
                                        z_scratch_chromosome,
                                        chrom,
                                        num_atoms,
                                        num_rotamers,
                                        fragments,
                                        frag_start_indices,
                                        frag_stop_indices);
      }
    }
  }

} // namespace mudock
