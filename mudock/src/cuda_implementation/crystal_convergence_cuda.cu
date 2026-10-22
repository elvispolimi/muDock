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
    const int ligand_id        = blockIdx.x;
    const int local_thread_id  = threadIdx.x;
    assert(blockDim.x == BLOCK_SIZE && warpSize == BLOCK_SIZE &&
           "Warpsize and the number of thread per block does not coincide");
    if (converged_ligands_b[ligand_id]) return;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // WARNING: CRITERION 1 AND 2 ARE ARTIFICAL CONVERGENCE CRITERION, NOT REAL ONE, BECAUSE THEY ASSUME
    //          WE HAVE THE CRYSTAL, WHICH IS OK FOR EVALUATING PERFORMANCE, BUT NOT USEFUL FOR REAL CASE
    //          SCENARIOS WHERE INSTEAD WE MUST CHECK CONVERGENCE LOOKING AT THE SCORE OR POPULATION
    //          STAGNATION ETC...
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // WARNING: TODO L: criterion on crystal_score is valid for only one ligand at the time
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // CRITERION 1: best_score <= (crystal_score + 1)
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // TODO L: implement 
    const fp_type *__restrict__ scores_l = scores_b + ligand_id * individuals_per_ligand;
    int best_index = 0;
    fp_type best_score = scores_l[best_index];
    for (int individual_index = 0; individual_index < individuals_per_ligand; individual_index++) {
      if (scores_l[individual_index] < best_score){
        best_index = individual_index;
        best_score = scores_l[individual_index];
      }
    }
    bool good_score = best_score <= (crystal_score + SCORE_TOLERANCE) ? true : false;


    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // CRITERION 2: RMSD best scoring pose < 1 Ångström
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    const int num_atoms   = num_atoms_b[ligand_id];
    const int atom_stride = ligand_id * batch_atoms;
    const fp_type* template_x_l = template_x_b + atom_stride;
    const fp_type* template_y_l = template_y_b + atom_stride;
    const fp_type* template_z_l = template_z_b + atom_stride;

    const fp_type *__restrict__ scratch_x = x_scratch_b + atom_stride * individuals_per_ligand;
    const fp_type *__restrict__ scratch_y = y_scratch_b + atom_stride * individuals_per_ligand;
    const fp_type *__restrict__ scratch_z = z_scratch_b + atom_stride * individuals_per_ligand;
    const fp_type *__restrict__ scratch_x_l = scratch_x + best_index * batch_atoms;
    const fp_type *__restrict__ scratch_y_l = scratch_y + best_index * batch_atoms;
    const fp_type *__restrict__ scratch_z_l = scratch_z + best_index * batch_atoms;

    fp_type sum_sq_dist = fp_type{0};
    for (int atom_index = 0; atom_index < num_atoms; ++atom_index) {
      const fp_type dx = template_x_l[atom_index] - scratch_x_l[atom_index];
      const fp_type dy = template_y_l[atom_index] - scratch_y_l[atom_index];
      const fp_type dz = template_z_l[atom_index] - scratch_z_l[atom_index];
      sum_sq_dist += dx * dx + dy * dy + dz * dz;
    }
    rmsd_best_pose_b[ligand_id] = std::sqrt(sum_sq_dist / static_cast<fp_type>(num_atoms));
    bool good_rmsd = rmsd_best_pose_b[ligand_id] < RMSD_TOLERANCE ? true : false;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // CRITERION 3: best score so far has not improved (significantly)
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    bool best_score_not_improved = (best_so_far_b[ligand_id] - best_score < best_score_diff_thld) ? true : false;
    
    if (best_score_not_improved) {
      for_how_long_best_b[ligand_id] += 1;
    } else {
      for_how_long_best_b[ligand_id] = 0;
    }
    
    bool stale_best_score = (for_how_long_best_b[ligand_id] >= tolerance_window) ? true : false;
    
    // Update best so far
    best_so_far_b[ligand_id] = std::min(best_so_far_b[ligand_id], best_score);

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // LOG progress
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // Mark convergence
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // TODO L: i moved autostop parameter here so that in every case we log rmsd, but technically not 
    //         needed. If autostop is inside genetic and decides if this convergence stage runs or not,
    //         it may not log it.
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    if (autostop && (good_score || good_rmsd || stale_best_score)) {
      converged_ligands_b[ligand_id] = generation;
    }

    // TODO L IMPORTANT: implement timer here as in cpp

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