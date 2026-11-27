#include <cassert>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/gh_implementation/geom_transform_gh.hpp>
#include <mudock/gh_implementation/mutate_gh.hpp>

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
        apply<cpu_vectorization::GH>(x_scratch_chromosome,
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
  void geometric<queue_gh>::operator()() {
    static constexpr char geom_region_name[] = "geometric_trasformation";
    this->get_queue()->invoke_kernel<geom_region_name>(geom_transform,
                                                       args.batch_ligands,
                                                       args.batch_atom_stride,
                                                       args.chromosomes_per_ligand,
                                                       args.chromosomes,
                                                       args.num_atoms,
                                                       args.num_rotamers,
                                                       args.original_x,
                                                       args.original_y,
                                                       args.original_z,
                                                       args.x_scratch,
                                                       args.y_scratch,
                                                       args.z_scratch,
                                                       args.fragments,
                                                       args.ligand_fragments_start,
                                                       args.frag_indices_start,
                                                       args.fragments_start_index,
                                                       args.fragments_stop_index);
  }

} // namespace mudock
