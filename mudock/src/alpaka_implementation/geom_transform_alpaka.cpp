#include <alpaka/alpaka.hpp>
#include <cmath>
#include <mudock/alpaka_implementation/geom_transform_alpaka.hpp>
#include <mudock/alpaka_implementation/invoke_kernel_alpaka.hpp>
#include <mudock/alpaka_implementation/mutate_alpaka.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/utils.hpp>
#ifndef MUDOCK_ALPAKA_BLOCK_SIZE
  #define MUDOCK_ALPAKA_BLOCK_SIZE 32
#endif

namespace mudock {
  namespace {
    template<int MAX_ATOMS>
    struct apply_alpaka {
      template<typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int chromosome_number,
                                    const int atom_stride,
                                    const fp_type* original_x,
                                    const fp_type* original_y,
                                    const fp_type* original_z,
                                    fp_type* scratch_x,
                                    fp_type* scratch_y,
                                    fp_type* scratch_z,
                                    const chromosome* chromosomes,
                                    const int* fragments,
                                    const int* ligand_fragments_start,
                                    const int* fragments_start_index,
                                    const int* fragments_stop_index,
                                    const int* frag_indices_start,
                                    const int* num_rotamers_b,
                                    const int* num_atoms_b) const {
        const int ligand_id = static_cast<int>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        const int thread_id = static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        const int thread_per_block =
            static_cast<int>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);

        const int num_atoms    = num_atoms_b[ligand_id];
        const int num_rotamers = num_rotamers_b[ligand_id];

        const fp_type* l_original_x         = original_x + ligand_id * atom_stride;
        const fp_type* l_original_y         = original_y + ligand_id * atom_stride;
        const fp_type* l_original_z         = original_z + ligand_id * atom_stride;
        fp_type* l_scratch_x                = scratch_x + ligand_id * atom_stride * chromosome_number;
        fp_type* l_scratch_y                = scratch_y + ligand_id * atom_stride * chromosome_number;
        fp_type* l_scratch_z                = scratch_z + ligand_id * atom_stride * chromosome_number;
        const chromosome* chromosomes_b     = chromosomes + ligand_id * chromosome_number;
        const auto* l_fragments             = fragments + ligand_fragments_start[ligand_id];
        const auto* l_frag_start_atom_index = fragments_start_index + frag_indices_start[ligand_id];
        const auto* l_frag_stop_atom_index  = fragments_stop_index + frag_indices_start[ligand_id];

        for (int chromosome_index = thread_id; chromosome_index < chromosome_number;
             chromosome_index += thread_per_block) {
          const chromosome& l_chromosomes = chromosomes_b[chromosome_index];
          fp_type* x_scratch_chromosome   = l_scratch_x + chromosome_index * atom_stride;
          fp_type* y_scratch_chromosome   = l_scratch_y + chromosome_index * atom_stride;
          fp_type* z_scratch_chromosome   = l_scratch_z + chromosome_index * atom_stride;

          for (int atom_index = 0; atom_index < MAX_ATOMS; ++atom_index) {
            if (atom_index < num_atoms) {
              x_scratch_chromosome[atom_index] = l_original_x[atom_index];
              y_scratch_chromosome[atom_index] = l_original_y[atom_index];
              z_scratch_chromosome[atom_index] = l_original_z[atom_index];
            }
          }

          translate_molecule_alpaka<MAX_ATOMS>(x_scratch_chromosome,
                                               y_scratch_chromosome,
                                               z_scratch_chromosome,
                                               l_chromosomes[0],
                                               l_chromosomes[1],
                                               l_chromosomes[2],
                                               num_atoms);
          rotate_molecule_alpaka<MAX_ATOMS>(x_scratch_chromosome,
                                            y_scratch_chromosome,
                                            z_scratch_chromosome,
                                            l_chromosomes[3],
                                            l_chromosomes[4],
                                            l_chromosomes[5],
                                            num_atoms);

          for (int i = 0; i < num_rotamers; ++i) {
            const int* bitmask = l_fragments + i * num_atoms;
            rotate_fragment_alpaka<MAX_ATOMS>(x_scratch_chromosome,
                                              y_scratch_chromosome,
                                              z_scratch_chromosome,
                                              bitmask,
                                              l_frag_start_atom_index[i],
                                              l_frag_stop_atom_index[i],
                                              l_chromosomes[6 + i],
                                              num_atoms);
          }
        }
      }
    };
  } // namespace

  template<>
  void geom_kernel<queue_alpaka>::operator()() {
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->invoke_kernel<apply_alpaka<max_atoms>>(batch_ligands,
                                                    MUDOCK_ALPAKA_BLOCK_SIZE,
                                                    chromsomes_per_ligand,
                                                    batch_atoms,
                                                    x_coords_b,
                                                    y_coords_b,
                                                    z_coords_b,
                                                    x_scratch_b,
                                                    y_scratch_b,
                                                    z_scratch_b,
                                                    chromosomes_b,
                                                    ligand_fragments_b,
                                                    ligand_fragments_start_b,
                                                    frag_start_indices_b,
                                                    frag_stop_indices_b,
                                                    frag_indices_start_b,
                                                    num_rotamers_b,
                                                    num_atoms_b);
        },
        batch_atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }
} // namespace mudock
