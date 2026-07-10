#include <cassert>
#include <mudock/cpp_implementation/geom_transform_cpp.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>

namespace mudock {
  namespace {
    void geom_transform(const int batch_ligands,
                        const int batch_atoms,
                        const int chromsomes_per_ligand,
                        const chromosome* __restrict__ chromosomes_b,
                        const int* __restrict__ num_atoms_b,
                        const int* __restrict__ num_rotamers_b,
                        const fp_type* __restrict__ x_coords_b,
                        const fp_type* __restrict__ y_coords_b,
                        const fp_type* __restrict__ z_coords_b,
                        fp_type* __restrict__ x_scratch_b,
                        fp_type* __restrict__ y_scratch_b,
                        fp_type* __restrict__ z_scratch_b,
                        const int* __restrict__ ligand_fragments_b,
                        const int* __restrict__ ligand_fragments_start_b,
                        const int* __restrict__ frag_indices_start_b,
                        const int* __restrict__ frag_start_indices_b,
                        const int* __restrict__ frag_stop_indices_b) {
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

      fp_type* __restrict__ x_scratch = x_scratch_b + stride_atoms * chromsomes_per_ligand;
      fp_type* __restrict__ y_scratch = y_scratch_b + stride_atoms * chromsomes_per_ligand;
      fp_type* __restrict__ z_scratch = z_scratch_b + stride_atoms * chromsomes_per_ligand;
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
        apply<cpu_vectorization::AUTO>(x_scratch_chromosome,
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
  } // namespace
  template<>
  void geom_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<geom_region_name>(geom_transform,
                                       batch_ligands,
                                       batch_atoms,
                                       chromsomes_per_ligand,
                                       chromosomes_b,
                                       num_atoms_b,
                                       num_rotamers_b,
                                       x_coords_b,
                                       y_coords_b,
                                       z_coords_b,
                                       x_scratch_b,
                                       y_scratch_b,
                                       z_scratch_b,
                                       ligand_fragments_b,
                                       ligand_fragments_start_b,
                                       frag_indices_start_b,
                                       frag_start_indices_b,
                                       frag_stop_indices_b);
  }
} // namespace mudock
