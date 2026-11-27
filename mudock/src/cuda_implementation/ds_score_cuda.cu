
#include <cuda_runtime.h>
#include <mudock/compute/devices_memory.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cuda_implementation/cuda_utils.cuh>
#include <mudock/cuda_implementation/ds_score_cuda.cuh>
#include <mudock/utils.hpp>

#define BUCKET_MULTIPLIER 2

namespace mudock {

  template<int MAX_ATOMS>
  __global__ void calc_drugscore(const int atom_stride,
                                 const int scores_per_ligand,
                                 const fp_type* __restrict__ batch_ligand_x,
                                 const fp_type* __restrict__ batch_ligand_y,
                                 const fp_type* __restrict__ batch_ligand_z,
                                 const ds_type* __restrict__ batch_ligand_ds_types,
                                 const int* batch_ligand_num_atoms,
                                 const fp_type* __restrict__ protein_x,
                                 const fp_type* __restrict__ protein_y,
                                 const fp_type* __restrict__ protein_z,
                                 const ds_type* __restrict__ protein_ds_types,
                                 const int protein_num_atoms,
                                 const fp_type* __restrict__ drugscore_scores,
                                 fp_type* __restrict__ batch_scores) {
    const unsigned int nd_idx     = ds2index(ds_type::ND);
    const unsigned int out_stride = outer_stride();
    const unsigned int int_stride = inner_stride();
    const unsigned int min_bkt    = min_bucket();
    const unsigned int max_bkt    = max_bucket();

    const int ligand_index    = blockIdx.x;
    const int local_thread_id = threadIdx.x;

    const int ligand_stride    = ligand_index * atom_stride * scores_per_ligand;
    const int ligand_num_atoms = batch_ligand_num_atoms[ligand_index];

    const fp_type* pool_ligand_x   = batch_ligand_x + ligand_stride;
    const fp_type* pool_ligand_y   = batch_ligand_y + ligand_stride;
    const fp_type* pool_ligand_z   = batch_ligand_z + ligand_stride;
    const ds_type* ligand_ds_types = batch_ligand_ds_types + ligand_index * atom_stride;

    fp_type* pool_scores = batch_scores + ligand_index * scores_per_ligand;

    for (int score_index{0}; score_index < scores_per_ligand; ++score_index) {
      const int pool_stride = score_index * atom_stride;

      const fp_type* ligand_x = pool_ligand_x + pool_stride;
      const fp_type* ligand_y = pool_ligand_y + pool_stride;
      const fp_type* ligand_z = pool_ligand_z + pool_stride;

      fp_type drugscore{0};

#pragma unroll
      for (int l_a_it = threadIdx.x; l_a_it < MAX_ATOMS; l_a_it += blockDim.x) {
        if (l_a_it < ligand_num_atoms) {
          // Ligand atom type
          const ds_type l_a_ds       = ligand_ds_types[l_a_it];
          const unsigned int l_a_idx = ds2index(l_a_ds);
          const bool l_a_valid       = l_a_idx != nd_idx;

          // Ligand atom coordinates
          const fp_type l_a_x = ligand_x[l_a_it];
          const fp_type l_a_y = ligand_y[l_a_it];
          const fp_type l_a_z = ligand_z[l_a_it];

          for (int p_a_it = 0; p_a_it < protein_num_atoms; ++p_a_it) {
            // Protein atom type
            const ds_type p_a_ds       = protein_ds_types[p_a_it];
            const unsigned int p_a_idx = ds2index(p_a_ds);
            const bool p_a_valid       = p_a_idx != nd_idx;

            // Protein atom coordinates
            const fp_type p_a_x = protein_x[p_a_it];
            const fp_type p_a_y = protein_y[p_a_it];
            const fp_type p_a_z = protein_z[p_a_it];

            // Distance calculation
            const fp_type dist_x = l_a_x - p_a_x;
            const fp_type dist_y = l_a_y - p_a_y;
            const fp_type dist_z = l_a_z - p_a_z;
            const fp_type d_two  = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;
            const fp_type dist   = sqrtf(d_two);

            // Bucketize
            const unsigned int dist_bucket = bucketize(dist);
            const bool range_valid         = (dist_bucket >= min_bkt) && (dist_bucket <= max_bkt);

            // Validity
            const bool valid = l_a_valid && p_a_valid && range_valid;

            // Get flat index
            const unsigned int bucket        = dist_bucket - min_bkt;
            const unsigned int temp_flat_idx = l_a_idx * out_stride + p_a_idx * int_stride + bucket;
            const unsigned int flat_idx      = temp_flat_idx * valid;

            // Get score
            const fp_type temp_score = drugscore_scores[flat_idx];
            const fp_type score      = temp_score * valid;

            drugscore += score;
          }
        }
      }
      // Intra block reduction (block is sized as a single warp for optimizations)
#pragma unroll
      for (int offset = BLOCK_SIZE / 2; offset > 0; offset /= 2) {
        drugscore += __shfl_down_sync(0xffffffff, drugscore, offset);
      }

      if (local_thread_id == 0) {
        pool_scores[score_index] = drugscore;
      }
    }
  }

  template<>
  void ds_score_kernel<queue_cuda>::operator()() {
    const int dev_id = q->get_id();

    void* args[] = {(void*) &batch_atoms,
                    (void*) &scores_per_ligand,
                    (void*) &batch_ligand_x,
                    (void*) &batch_ligand_y,
                    (void*) &batch_ligand_z,
                    (void*) &batch_ligand_ds_types,
                    (void*) &batch_ligand_num_atoms,
                    (void*) &protein_x,
                    (void*) &protein_y,
                    (void*) &protein_z,
                    (void*) &protein_ds_types,
                    (void*) &protein_num_atoms,
                    (void*) &ds_scores,
                    (void*) &scores_b};

    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->launch_kernel((void*) calc_drugscore<max_atoms>, args, batch_ligands, BLOCK_SIZE);
        },
        batch_atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

  //===------------------------------------------------------------------------------------------------------
  // PRECOMPUTED GRID-BASED KERNEL
  //===------------------------------------------------------------------------------------------------------
  __device__ __constant__ fp_type ds_map_min_const[3];
  __device__ __constant__ fp_type ds_map_max_const[3];

  constexpr int ds_k_max_devices = 16;
  device_memory_array<ds_k_max_devices, fp_type> ds_cuda_constant_memory;

  void init_device(const int dev,
                   const fp_type* map_min,
                   const fp_type* map_max) {
    // Thread-safe, exactly-once init per device:
    ds_cuda_constant_memory.init(dev, std::function<std::unique_ptr<fp_type>()>([&]() {
                                   MUDOCK_CHECK(cudaSetDevice(dev));
                                   // MUDOCK_CHECK(cudaStreamCreate(&stream));
                                   MUDOCK_CHECK(cudaMemcpyToSymbol(ds_map_min_const,
                                                                   map_min,
                                                                   3 * sizeof(fp_type),
                                                                   0,
                                                                   cudaMemcpyDeviceToDevice));
                                   MUDOCK_CHECK(cudaMemcpyToSymbol(ds_map_max_const,
                                                                   map_max,
                                                                   3 * sizeof(fp_type),
                                                                   0,
                                                                   cudaMemcpyDeviceToDevice));
                                   MUDOCK_CHECK(cudaDeviceSynchronize());
                                   // TODO implement RAI to release constant memory
                                   return std::make_unique<fp_type>(0);
                                 }));
  }

  template<int MAX_ATOMS>
  __global__ void calc_drugscore_precomp(const int atom_stride,
                                         const int scores_per_ligand,
                                         const fp_type* __restrict__ batch_ligand_x,
                                         const fp_type* __restrict__ batch_ligand_y,
                                         const fp_type* __restrict__ batch_ligand_z,
                                         const ds_type* __restrict__ batch_ligand_ds_types,
                                         const int* batch_ligand_num_atoms,
                                         const int map_index_x,
                                         const int map_index_xy,
                                         const int map_index_xyz,
                                         const fp_type* __restrict__ ds_grids,
                                         fp_type* __restrict__ batch_scores) {
    const unsigned int nd_idx = ds2index(ds_type::ND);

    const int ligand_index    = blockIdx.x;
    const int local_thread_id = threadIdx.x;

    const int ligand_stride    = ligand_index * atom_stride * scores_per_ligand;
    const int ligand_num_atoms = batch_ligand_num_atoms[ligand_index];

    const fp_type* pool_ligand_x   = batch_ligand_x + ligand_stride;
    const fp_type* pool_ligand_y   = batch_ligand_y + ligand_stride;
    const fp_type* pool_ligand_z   = batch_ligand_z + ligand_stride;
    const ds_type* ligand_ds_types = batch_ligand_ds_types + ligand_index * atom_stride;

    fp_type* pool_scores = batch_scores + ligand_index * scores_per_ligand;

    for (int score_index{0}; score_index < scores_per_ligand; ++score_index) {
      const int pool_stride = score_index * atom_stride;

      const fp_type* ligand_x = pool_ligand_x + pool_stride;
      const fp_type* ligand_y = pool_ligand_y + pool_stride;
      const fp_type* ligand_z = pool_ligand_z + pool_stride;

      fp_type drugscore{0};

#pragma unroll
      for (int l_a_it = threadIdx.x; l_a_it < MAX_ATOMS; l_a_it += blockDim.x) {
        if (l_a_it < ligand_num_atoms) {
          // Ligand atom type
          const ds_type l_a_ds       = ligand_ds_types[l_a_it];
          const unsigned int l_a_idx = ds2index(l_a_ds);
          const bool l_a_valid       = l_a_idx != nd_idx;

          // Ligand atom coordinates
          const fp_type l_a_x = ligand_x[l_a_it];
          const fp_type l_a_y = ligand_y[l_a_it];
          const fp_type l_a_z = ligand_z[l_a_it];

          // Compute bounds check
          const bool in_bounds_x  = (l_a_x >= ds_map_min_const[0]) && (l_a_x <= ds_map_max_const[0]);
          const bool in_bounds_y  = (l_a_y >= ds_map_min_const[1]) && (l_a_y <= ds_map_max_const[1]);
          const bool in_bounds_z  = (l_a_z >= ds_map_min_const[2]) && (l_a_z <= ds_map_max_const[2]);
          const bool valid_bounds = in_bounds_x && in_bounds_y && in_bounds_z;

          // Compute texture coordinates
          const fp_type grid_coord_x = (l_a_x - ds_map_min_const[0]) * inv_spacing;
          const fp_type grid_coord_y = (l_a_y - ds_map_min_const[1]) * inv_spacing;
          const fp_type grid_coord_z = (l_a_z - ds_map_min_const[2]) * inv_spacing;

          const int atom_index_x  = static_cast<int>(grid_coord_x + fp_type(0.5));
          const fp_type weight_x0 = grid_coord_x - static_cast<fp_type>(atom_index_x);
          const fp_type weight_x1 = fp_type(1) - weight_x0;

          const int atom_index_y  = static_cast<int>(grid_coord_y + fp_type(0.5));
          const fp_type weight_y0 = grid_coord_y - static_cast<fp_type>(atom_index_y);
          const fp_type weight_y1 = fp_type(1) - weight_y0;

          const int atom_index_z  = static_cast<int>(grid_coord_z + fp_type(0.5));
          const fp_type weight_z0 = grid_coord_z - static_cast<fp_type>(atom_index_z);
          const fp_type weight_z1 = fp_type(1) - weight_z0;

          // Compute coefficients
          const fp_type c0 = weight_x1 * weight_y1 * weight_z1;
          const fp_type c1 = weight_x1 * weight_y1 * weight_z0;
          const fp_type c2 = weight_x1 * weight_y0 * weight_z1;
          const fp_type c3 = weight_x1 * weight_y0 * weight_z0;
          const fp_type c4 = weight_x0 * weight_y1 * weight_z1;
          const fp_type c5 = weight_x0 * weight_y1 * weight_z0;
          const fp_type c6 = weight_x0 * weight_y0 * weight_z1;
          const fp_type c7 = weight_x0 * weight_y0 * weight_z0;

          const int base_index = map_index_xy * atom_index_z + atom_index_y * map_index_x + atom_index_x;
          const int map_offset = (base_index + map_index_xyz * l_a_idx) * l_a_valid; // safe-zero if not valid

          // Inline trilinear interpolation for better vectorization
          fp_type value = fp_type(0);
          value += c0 * ds_grids[map_offset];
          value += c1 * ds_grids[map_offset + map_index_xy];
          value += c2 * ds_grids[map_offset + map_index_x];
          value += c3 * ds_grids[map_offset + map_index_x + map_index_xy];
          value += c4 * ds_grids[map_offset + 1];
          value += c5 * ds_grids[map_offset + 1 + map_index_xy];
          value += c6 * ds_grids[map_offset + 1 + map_index_x];
          value += c7 * ds_grids[map_offset + 1 + map_index_x + map_index_xy];

          drugscore += value * (l_a_valid & valid_bounds);
        }
      }
      // Intra block reduction (block is sized as a single warp for optimizations)
#pragma unroll
      for (int offset = BLOCK_SIZE / 2; offset > 0; offset /= 2) {
        drugscore += __shfl_down_sync(0xffffffff, drugscore, offset);
      }

      if (local_thread_id == 0) {
        pool_scores[score_index] = drugscore;
      }
    }
  }

  template<>
  void ds_score_precomp_kernel<queue_cuda>::operator()() {
    const int dev_id = q->get_id();
    init_device(dev_id, minimum, maximum);

    void* args[] = {(void*) &batch_atoms,
                    (void*) &scores_per_ligand,
                    (void*) &batch_ligand_x,
                    (void*) &batch_ligand_y,
                    (void*) &batch_ligand_z,
                    (void*) &batch_ligand_ds_types,
                    (void*) &batch_ligand_num_atoms,
                    (void*) &map_index_x,
                    (void*) &map_index_xy,
                    (void*) &map_index_xyz,                    
                    (void*) &grid_maps,
                    (void*) &scores_b};

    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->launch_kernel((void*) calc_drugscore_precomp<max_atoms>, args, batch_ligands, BLOCK_SIZE);
        },
        batch_atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

  //===------------------------------------------------------------------------------------------------------
  // BATCH SIZE COMPUTATION
  //===------------------------------------------------------------------------------------------------------
  template<int MAX_ATOMS>
  batch_multiple get_ds_evaluate_fitness_batch(const int device_id) {
    return get_kernel_batch_multiple_cuda<calc_drugscore<MAX_ATOMS>>(device_id,
                                                                  BLOCK_SIZE,
                                                                  0,
                                                                  "adt_score::calc_energy");
  }

  template<>
  batch_multiple get_ds_score_batch_multiple<queue_cuda>(const int atoms, std::shared_ptr<queue_cuda> q_b) {
    batch_multiple bucket_multiple{};
    const int device_id = q_b->get_id();
    constexpr_for<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>([&](const auto atom_index) {
      const auto n_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
      if (atoms == n_atoms)
        bucket_multiple = get_ds_evaluate_fitness_batch<n_atoms>(device_id);
    });
    if (bucket_multiple.total_multiple() <= 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");
    return normalize_batch_multiple(bucket_multiple);
  };

}; // namespace mudock
