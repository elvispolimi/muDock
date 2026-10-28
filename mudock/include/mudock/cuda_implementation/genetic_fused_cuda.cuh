#pragma once

#include <mudock/compute/batch_multiple.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/type_alias.hpp>
#include <stdexcept>
#include <memory>
#include <cstddef>
#ifdef MUDOCK_USE_CUDA
  #include <mudock/cuda_implementation/queue_cuda.cuh>
#endif

// Source-level experiment switches for the fused genetic kernel.
// Uncomment to process one ligand per block with one warp instead of four.
// #define FUSED_GENETIC_SINGLE_WARP
#ifdef FUSED_GENETIC_SINGLE_WARP
  #define FUSED_WARPS_PER_BLOCK 1
#else
  #define FUSED_WARPS_PER_BLOCK 4
#endif

// Comment out to store per-warp pose coordinates in the global scratch slab.
#define FUSED_GENETIC_USE_SHARED_MEMORY

// Comment out to use a runtime rotamer loop instead of template specialization.
#define FUSED_GENETIC_TEMPLATE_ROTAMERS
// 1 requests full unrolling; larger values reduce the template unroll factor.
#define FUSED_ROTAMER_UNROLL_DIVISOR 1

namespace mudock {
  struct genetic_fused_data {
    int ligands;
    int atoms;
    int max_rotamers;
    int population;
    int generations;
    int tournament_length;
    fp_type mutation_probability;
    std::size_t seed;
    fp_type* current;
    fp_type* next;
    fp_type* scores;
    fp_type* best_scores;
    fp_type* best_chromosomes;
    fp_type* pose_scratch;
    const int* num_atoms;
    const int* num_rotamers;
    const fp_type* x_coords;
    const fp_type* y_coords;
    const fp_type* z_coords;
    const int* fragments;
    const int* fragments_start;
    const int* frag_indices_start;
    const int* frag_start;
    const int* frag_stop;
    const int* num_nonbond;
    const fp_type* vols;
    const fp_type* solpars;
    const fp_type* charges;
    const int* map_offsets;
    const int* nonbond_a1;
    const int* nonbond_a2;
    const fp_type* nonbond_cA;
    const fp_type* nonbond_cB;
    const int* nonbond_xB;
    const fp_type* grid_maps;
    const fp_type* map_min;
    const fp_type* map_max;
    const fp_type* map_center;
    int map_index_x;
    int map_index_xy;
    int map_index_xyz;
  };

  template<typename queue_t>
  struct genetic_fused_launcher {
    genetic_fused_launcher(genetic_fused_data data_, std::shared_ptr<queue_t>) : data(data_) {}
    void operator()() { throw std::runtime_error("genetic_fused is supported only by CUDA"); }
  private:
    genetic_fused_data data;
  };

#ifdef MUDOCK_USE_CUDA
  batch_multiple get_genetic_fused_batch_multiple_cuda(int device_id, int atoms, int max_rotamers);

  template<>
  struct genetic_fused_launcher<queue_cuda> {
    genetic_fused_launcher(genetic_fused_data data_, std::shared_ptr<queue_cuda> queue_)
        : data(data_), q(std::move(queue_)) {}
    void operator()();
  private:
    genetic_fused_data data;
    std::shared_ptr<queue_cuda> q;
  };
#endif
} // namespace mudock
