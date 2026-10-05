/**
 * @file geom_transform_alpaka.cpp
 * @brief Implementation of the geometric transformation kernel for Alpaka.
 * @details Implements the `apply_alpaka` kernel and batch multiple calculation
 *          for evaluating ligand conformations by applying genotype translation,
 *          quaternion/Euler rotations, and flexible fragment torsions.
 */

#include <alpaka/alpaka.hpp>
#include <cmath>
#include <mudock/alpaka_implementation/geom_transform_alpaka.hpp>
#include <mudock/alpaka_implementation/invoke_kernel_alpaka.hpp>
#include <mudock/alpaka_implementation/mutate_alpaka.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/utils.hpp>


namespace mudock {
  namespace {
    /**
     * @struct apply_alpaka
     * @brief Device kernel performing geometric conformation transformations on molecular candidates.
     * @details For each candidate ligand chromosome, reads atomic baseline coordinates, applies
     *          global 3D translation and centroid rotation, followed by internal rotatable bond
     *          torsions on molecular sub-fragments according to bitmasks.
     *
     * @tparam MAX_ATOMS Static upper bound on atom count used for unrolling inner loops.
     */
    template<int MAX_ATOMS>
    struct apply_alpaka {
      /**
       * @brief Kernel body executing translation, rotation, and fragment torsion transformations.
       * @tparam TAcc Alpaka accelerator type.
       * @param[in] acc Reference to the execution context.
       * @param chromosome_number Number of candidate chromosomes per ligand pose.
       * @param atom_stride Stride in elements separating coordinate planes across molecules.
       * @param[in] original_x Baseline X atomic coordinates.
       * @param[in] original_y Baseline Y atomic coordinates.
       * @param[in] original_z Baseline Z atomic coordinates.
       * @param[out] scratch_x Scratchpad destination buffer for transformed X coordinates.
       * @param[out] scratch_y Scratchpad destination buffer for transformed Y coordinates.
       * @param[out] scratch_z Scratchpad destination buffer for transformed Z coordinates.
       * @param[in] chromosomes Candidate chromosomes containing rotation and torsion angles.
       * @param[in] fragments Array of fragment membership bitmasks.
       * @param[in] ligand_fragments_start Array of starting offsets in fragments array per ligand.
       * @param[in] fragments_start_index Indices of origin atoms for fragment rotation axes.
       * @param[in] fragments_stop_index Indices of terminal atoms for fragment rotation axes.
       * @param[in] frag_indices_start Array of offsets into fragment axis index arrays per ligand.
       * @param[in] num_rotamers_b Number of rotatable bonds per ligand.
       * @param[in] num_atoms_b Number of atoms per ligand.
       */
      template<typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int chromosome_number,
                                    const int atom_stride,
                                    const fp_type* __restrict__ original_x,
                                    const fp_type* __restrict__ original_y,
                                    const fp_type* __restrict__ original_z,
                                    fp_type* __restrict__ scratch_x,
                                    fp_type* __restrict__ scratch_y,
                                    fp_type* __restrict__ scratch_z,
                                    const chromosome* __restrict__ chromosomes,
                                    const int* __restrict__ fragments,
                                    const int* __restrict__ ligand_fragments_start,
                                    const int* __restrict__ fragments_start_index,
                                    const int* __restrict__ fragments_stop_index,
                                    const int* __restrict__ frag_indices_start,
                                    const int* __restrict__ num_rotamers_b,
                                    const int* __restrict__ num_atoms_b) const {
        const int ligand_id = static_cast<int>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        const int thread_id = static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);

        const int num_atoms    = num_atoms_b[ligand_id];
        const int num_rotamers = num_rotamers_b[ligand_id];

        const fp_type* __restrict__ l_original_x    = original_x + ligand_id * atom_stride;
        const fp_type* __restrict__ l_original_y    = original_y + ligand_id * atom_stride;
        const fp_type* __restrict__ l_original_z    = original_z + ligand_id * atom_stride;
        fp_type* __restrict__ l_scratch_x           = scratch_x + ligand_id * atom_stride * chromosome_number;
        fp_type* __restrict__ l_scratch_y           = scratch_y + ligand_id * atom_stride * chromosome_number;
        fp_type* __restrict__ l_scratch_z           = scratch_z + ligand_id * atom_stride * chromosome_number;
        const chromosome* chromosomes_b             = chromosomes + ligand_id * chromosome_number;
        const auto* __restrict__ l_fragments        = fragments + ligand_fragments_start[ligand_id];
        const auto* __restrict__ l_frag_start_atom_index = fragments_start_index + frag_indices_start[ligand_id];
        const auto* __restrict__ l_frag_stop_atom_index  = fragments_stop_index + frag_indices_start[ligand_id];

        for (int chromosome_index = 0; chromosome_index < chromosome_number; ++chromosome_index) {
          const chromosome& l_chromosomes            = chromosomes_b[chromosome_index];
          fp_type* __restrict__ x_scratch_chromosome = l_scratch_x + chromosome_index * atom_stride;
          fp_type* __restrict__ y_scratch_chromosome = l_scratch_y + chromosome_index * atom_stride;
          fp_type* __restrict__ z_scratch_chromosome = l_scratch_z + chromosome_index * atom_stride;

          ALPAKA_UNROLL(MUDOCK_ATOM_LOOP_UNROLL_FACTOR(MAX_ATOMS, MUDOCK_ALPAKA_BLOCK_SIZE))
          for (int i = 0; i < MAX_ATOMS; i += MUDOCK_ALPAKA_BLOCK_SIZE) {
            const int atom_index = i + thread_id;
            if (atom_index < num_atoms) {
              x_scratch_chromosome[atom_index] = l_original_x[atom_index];
              y_scratch_chromosome[atom_index] = l_original_y[atom_index];
              z_scratch_chromosome[atom_index] = l_original_z[atom_index];
            }
          }

          translate_molecule_alpaka<MAX_ATOMS, MUDOCK_ALPAKA_BLOCK_SIZE>(acc,
                                                                         x_scratch_chromosome,
                                                                         y_scratch_chromosome,
                                                                         z_scratch_chromosome,
                                                                         l_chromosomes[0],
                                                                         l_chromosomes[1],
                                                                         l_chromosomes[2],
                                                                         num_atoms);

          rotate_molecule_alpaka<MAX_ATOMS, MUDOCK_ALPAKA_BLOCK_SIZE>(acc,
                                                                      x_scratch_chromosome,
                                                                      y_scratch_chromosome,
                                                                      z_scratch_chromosome,
                                                                      l_chromosomes[3],
                                                                      l_chromosomes[4],
                                                                      l_chromosomes[5],
                                                                      num_atoms);
          alpaka::syncBlockThreads(acc);

          ALPAKA_UNROLL(MUDOCK_UNROLL_FACTOR)
          for (int i = 0; i < num_rotamers; ++i) {
            const int* __restrict__ bitmask = l_fragments + i * num_atoms;
            rotate_fragment_alpaka<MAX_ATOMS, MUDOCK_ALPAKA_BLOCK_SIZE>(acc,
                                                                        x_scratch_chromosome,
                                                                        y_scratch_chromosome,
                                                                        z_scratch_chromosome,
                                                                        bitmask,
                                                                        l_frag_start_atom_index[i],
                                                                        l_frag_stop_atom_index[i],
                                                                        l_chromosomes[6 + i],
                                                                        num_atoms);
            alpaka::syncBlockThreads(acc);
          }
        }
      }
    };
  } // namespace

  /**
   * @brief Dispatches the geometric transformation kernel for the current molecular batch.
   * @details Switches at compile-time between static atom count clusters (`MAX_ATOMS`) to
   *          instantiate the optimal unrolled template instance of `apply_alpaka`.
   */
  template<>
  void geom_kernel<queue_alpaka>::operator()() {
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->invoke_kernel<apply_alpaka<max_atoms>>(batch_ligands,
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

  /**
   * @brief Calculates hardware-aware batch multiple alignment for geometric transformation.
   * @details Inspects the accelerator device's multiprocessor count (`num_sms`) and heuristics
   *          based on atom count to determine optimal block occupancy per streaming multiprocessor.
   * @param atoms Target atom cluster size.
   * @param q_b Pointer to the Alpaka execution queue.
   * @return Normalized batch_multiple structure specifying scheduling constraints.
   */
  template<>
  batch_multiple get_geom_transform_batch_multiple<queue_alpaka>(const int atoms,
                                                                 std::shared_ptr<queue_alpaka> q_b) {
    batch_multiple bucket_multiple{};
    const auto& dev = q_b->native_device();
    const int num_sms = static_cast<int>(alpaka::getAccDevProps<alpaka_backend::acc>(dev).m_multiProcessorCount);

    constexpr_for<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto n_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          if (atoms == n_atoms) {
            if constexpr (MUDOCK_ALPAKA_BLOCK_SIZE == 1) {
              // CPU backends
              bucket_multiple = {1, num_sms};
            } else {
              // GPU backends
              int blocks_per_sm = 16;
              if constexpr (n_atoms > 64)  blocks_per_sm = 12;
              if constexpr (n_atoms > 128) blocks_per_sm = 8;
              if constexpr (n_atoms > 192) blocks_per_sm = 4;
              bucket_multiple = {blocks_per_sm, num_sms};
            }
          }
        });

    if (bucket_multiple.total_multiple() <= 0)
      throw std::runtime_error("Compilation error: there is a bucket of atoms number which it is not handled.");

    mudock::info("ALPAKA GEOM batch multiple for ", atoms, " atoms -> ", bucket_multiple.total_multiple());
    return normalize_batch_multiple(bucket_multiple);
  }
} // namespace mudock
