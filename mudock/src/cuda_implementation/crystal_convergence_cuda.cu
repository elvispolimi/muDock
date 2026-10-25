#include <mudock/type_alias.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/molecule/constraints.hpp>
#include <mudock/cuda_implementation/queue_cuda.cuh>
#include <mudock/compute/crystal_convergence.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cuda_implementation/crystal_convergence_cuda.cuh>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/molecule.hpp>
#include <mudock/utils.hpp>
#include <mudock/cuda_implementation/cuda_utils.cuh>
#include <cmath>

namespace mudock {

  __device__ __forceinline__ fp_type warp_allreduce_sum(fp_type value) {
    MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
    for (int offset = BLOCK_SIZE / 2; offset > 0; offset /= 2) {
      value += __shfl_xor_sync(0xffffffff, value, offset);
    }
    return value;
  }

  template<int MAX_ATOMS>
  __global__ void check_convergence_gpu(const int batch_ligands,
                                        const int batch_atoms,
                                        const bool autostop,
                                        const int individuals_per_ligand,
                                        const int* __restrict__ num_atoms_b,
                                        const int generation,
                                        int* __restrict__ converged_ligands_b,
                                        const fp_type* __restrict__ template_x_b,
                                        const fp_type* __restrict__ template_y_b,
                                        const fp_type* __restrict__ template_z_b,
                                        const fp_type *__restrict__ x_scratch_b,
                                        const fp_type *__restrict__ y_scratch_b,
                                        const fp_type *__restrict__ z_scratch_b,
                                        fp_type *__restrict__ rmsd_best_pose_b,
                                        const fp_type *__restrict__ scores_b,
                                        const fp_type crystal_score,
                                        fp_type* __restrict__ best_so_far_b,
                                        int* __restrict__ for_how_long_best_b,
                                        const int tolerance_window,
                                        const fp_type best_score_diff_thld) {
    const int ligand_id       = blockIdx.x;
    const int local_thread_id = threadIdx.x;
    assert(blockDim.x == BLOCK_SIZE && warpSize == BLOCK_SIZE &&
           "Warpsize and the number of thread per block does not coincide");

    // Block-uniform early exit: all lanes leave together, so the full-mask shuffles below are safe.
    if (converged_ligands_b[ligand_id]) return;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // WARNING: CRITERION 1 AND 2 ARE ARTIFICAL CONVERGENCE CRITERION, NOT REAL ONE, BECAUSE THEY ASSUME
    //          WE HAVE THE CRYSTAL, WHICH IS OK FOR EVALUATING PERFORMANCE, BUT NOT USEFUL FOR REAL CASE
    //          SCENARIOS WHERE INSTEAD WE MUST CHECK CONVERGENCE LOOKING AT THE SCORE OR POPULATION
    //          STAGNATION ETC...
    // WARNING: TODO L: criterion on crystal_score is valid for only one ligand at the time
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    constexpr int no_index = std::numeric_limits<int>::max();
    const fp_type *__restrict__ scores_l = scores_b + ligand_id * individuals_per_ligand;

    int best_index     = no_index;
    fp_type best_score = std::numeric_limits<fp_type>::infinity();
    for (int i = local_thread_id; i < individuals_per_ligand; i += BLOCK_SIZE) {
      const fp_type s = scores_l[i];
      if (best_index == no_index || s < best_score) { // i increases per thread: first minimum is kept
        best_score = s;
        best_index = i;
      }
    }
    MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
    for (int offset = BLOCK_SIZE / 2; offset > 0; offset /= 2) {
      const fp_type other_score = __shfl_xor_sync(0xFFFFFFFF, best_score, offset);
      const int other_index     = __shfl_xor_sync(0xFFFFFFFF, best_index, offset);
      if (other_index != no_index &&
          (best_index == no_index || other_score < best_score ||
           (other_score == best_score && other_index < best_index))) {
        best_score = other_score;
        best_index = other_index;
      }
    }
    // xor-butterfly: every lane now holds the same (best_score, best_index).

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // CRITERION 2: RMSD of the best scoring pose, atoms scattered over the threads
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    const int num_atoms         = num_atoms_b[ligand_id];
    const int atom_stride       = ligand_id * batch_atoms;
    const fp_type* template_x_l = template_x_b + atom_stride;
    const fp_type* template_y_l = template_y_b + atom_stride;
    const fp_type* template_z_l = template_z_b + atom_stride;

    const fp_type *__restrict__ scratch_x_l = x_scratch_b + atom_stride * individuals_per_ligand + best_index * batch_atoms;
    const fp_type *__restrict__ scratch_y_l = y_scratch_b + atom_stride * individuals_per_ligand + best_index * batch_atoms;
    const fp_type *__restrict__ scratch_z_l = z_scratch_b + atom_stride * individuals_per_ligand + best_index * batch_atoms;

    fp_type sum_sq_dist = fp_type{0};
    for (int atom_index = local_thread_id; atom_index < num_atoms; atom_index += BLOCK_SIZE) {
      const fp_type dx = template_x_l[atom_index] - scratch_x_l[atom_index];
      const fp_type dy = template_y_l[atom_index] - scratch_y_l[atom_index];
      const fp_type dz = template_z_l[atom_index] - scratch_z_l[atom_index];
      sum_sq_dist += dx * dx + dy * dy + dz * dz;
    }
    sum_sq_dist = warp_allreduce_sum(sum_sq_dist); // all lanes call it (no divergence before)

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // Everything below reads/modifies the per-ligand persistent state: a SINGLE lane does it
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    if (local_thread_id != 0) return;

    const fp_type rmsd = std::sqrt(sum_sq_dist / static_cast<fp_type>(num_atoms));
    rmsd_best_pose_b[ligand_id] = rmsd;

    // CRITERION 1: best_score <= crystal_score + tolerance
    const bool good_score = best_score <= (crystal_score + SCORE_TOLERANCE);
    // CRITERION 2: RMSD < tolerance
    const bool good_rmsd = rmsd < RMSD_TOLERANCE;

    // CRITERION 3: best score so far has not improved (significantly)
    const fp_type best_so_far = best_so_far_b[ligand_id];
    int stale_count           = for_how_long_best_b[ligand_id];
    stale_count               = (best_so_far - best_score < best_score_diff_thld) ? stale_count + 1 : 0;
    for_how_long_best_b[ligand_id] = stale_count;
    const bool stale_best_score    = stale_count >= tolerance_window;

    best_so_far_b[ligand_id] = std::min(best_so_far, best_score);

    // Mark convergence
    if (autostop && (good_score || good_rmsd || stale_best_score)) {
      converged_ligands_b[ligand_id] = generation;
    }

    // TODO L IMPORTANT: implement timer here as in cpp, if needed

  }


  template<>
  void crystal_convergence_kernel<queue_cuda>::operator()() {
    void* args[] = {(void*) &batch_ligands,
                    (void*) &batch_atoms,
                    (void*) &autostop,
                    (void*) &population_number,
                    (void*) &num_atoms_b,
                    (void*) &current_generation,
                    (void*) &converged_ligands_b,
                    (void*) &template_x_b,
                    (void*) &template_y_b,
                    (void*) &template_z_b,
                    (void*) &x_scratch_b,
                    (void*) &y_scratch_b,
                    (void*) &z_scratch_b,
                    (void*) &rmsd_best_pose_b,
                    (void*) &scores_b,
                    (void*) &crystal_score,
                    (void*) &best_so_far_b,
                    (void*) &for_how_long_best_b,
                    (void*) &tolerance_window,
                    (void*) &best_score_diff_thld};
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
      [&](const auto atom_index) {
        const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
        q->launch_kernel((void*) check_convergence_gpu<max_atoms>, args, batch_ligands, BLOCK_SIZE);
      },
      batch_atoms,
      reorder_buffer<static_molecule>::atoms_clusters.data());
    current_generation++;

  }

  template<int MAX_ATOMS>
  batch_multiple get_convergence_batch(const int device_id) {
    return get_kernel_batch_multiple_cuda<check_convergence_gpu<MAX_ATOMS>>(device_id,
                                                                            BLOCK_SIZE,
                                                                            0,
                                                                            "convergence::check_convergence_gpu");
  }

  template<>
  batch_multiple get_crystal_conv_batch_multiple<queue_cuda>(const int atoms, std::shared_ptr<queue_cuda> q_b) {
    batch_multiple bucket_multiple{};
    const int device_id = q_b->get_id();
    constexpr_for<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>([&](const auto atom_index) {
      const auto n_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
      if (atoms == n_atoms)
        bucket_multiple = get_convergence_batch<n_atoms>(device_id);
    });
    if (bucket_multiple.total_multiple() <= 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");
    return normalize_batch_multiple(bucket_multiple);
  }

} // namespace mudock