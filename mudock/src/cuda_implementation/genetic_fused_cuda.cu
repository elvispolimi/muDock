#include <cmath>
#include <limits>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cuda_implementation/adt_score_cuda.cuh>
#include <mudock/cuda_implementation/cuda_random.cuh>
#include <mudock/cuda_implementation/cuda_utils.cuh>
#include <mudock/cuda_implementation/genetic_fused_cuda.cuh>
#include <mudock/cuda_implementation/mutate.cuh>
#include <mudock/utils.hpp>

#define FUSED_FLATTENED_3D(x, y, z, ix, ixy) ((ixy) * (z) + (y) * (ix) + (x))
#define FUSED_WARPS_PER_BLOCK                4
#define FUSED_BLOCK_SIZE                     (BLOCK_SIZE * FUSED_WARPS_PER_BLOCK)
// 1 requests full unrolling; larger values reduce the template unroll factor.
#define FUSED_ROTAMER_UNROLL_DIVISOR 1
#if defined(MUDOCK_GENETIC_FUSED_LAUNCH_BOUNDS_ATOM_SCALED_MAX)
  #define FUSED_GENETIC_LAUNCH_BOUNDS(max_atoms) \
    __launch_bounds__(FUSED_BLOCK_SIZE, (max_atoms) >= 160 ? 1 : ((max_atoms) >= 128 ? 2 : 8))
#elif defined(MUDOCK_GENETIC_FUSED_LAUNCH_BOUNDS_ATOM_SCALED)
  #define FUSED_GENETIC_LAUNCH_BOUNDS(max_atoms) \
    __launch_bounds__(FUSED_BLOCK_SIZE, (max_atoms) >= 160 ? 2 : ((max_atoms) >= 128 ? 4 : 8))
#else
  #define FUSED_GENETIC_LAUNCH_BOUNDS(max_atoms)
#endif

namespace mudock {
  __device__ __forceinline__ fp_type fused_random(curandState& state, const fp_type min, const fp_type max) {
    const fp_type value = is_debug() ? fp_type{0.4} : curand_uniform(&state);
    return static_cast<fp_type>(value * (max - min) + min);
  }

  __device__ __forceinline__ int fused_random_int(curandState& state, const int min, const int max) {
    return static_cast<int>(fused_random(state, static_cast<fp_type>(min), static_cast<fp_type>(max)));
  }

  __device__ __forceinline__ int fused_select(curandState& state,
                                              const int population,
                                              const int tournament_length,
                                              const fp_type* __restrict__ scores) {
    int best = fused_random_int(state, 0, population - 1);
    for (int i = 0; i < tournament_length; ++i) {
      const int candidate = fused_random_int(state, 0, population - 1);
      if (scores[candidate] < scores[best])
        best = candidate;
    }
    return best;
  }

  __device__ __forceinline__ fp_type fused_trilinear(const fp_type* __restrict__ map,
                                                     const fp_type* __restrict__ coeffs,
                                                     const int ix,
                                                     const int ixy) {
    fp_type value = 0;
    value         = coeffs[0] * map[0] + value;
    value         = coeffs[1] * map[ixy] + value;
    value         = coeffs[2] * map[ix] + value;
    value         = coeffs[3] * map[ix + ixy] + value;
    value         = coeffs[4] * map[1] + value;
    value         = coeffs[5] * map[1 + ixy] + value;
    value         = coeffs[6] * map[1 + ix] + value;
    value         = coeffs[7] * map[1 + ix + ixy] + value;
    return value;
  }

  template<int MAX_ATOMS>
  __device__ __forceinline__ fp_type fused_score_pose(const int ligand,
                                                      const int lane,
                                                      const int atoms,
                                                      const int* __restrict__ atom_counts,
                                                      const int* __restrict__ rotamer_counts,
                                                      const fp_type* __restrict__ x_scratch,
                                                      const fp_type* __restrict__ y_scratch,
                                                      const fp_type* __restrict__ z_scratch,
                                                      const fp_type* __restrict__ map_min,
                                                      const fp_type* __restrict__ map_max,
                                                      const fp_type* __restrict__ map_center,
                                                      const fp_type* __restrict__ vols,
                                                      const fp_type* __restrict__ solpars,
                                                      const fp_type* __restrict__ charges,
                                                      const int* __restrict__ map_offsets,
                                                      const int* __restrict__ num_nonbond,
                                                      const int* __restrict__ nonbond_a1,
                                                      const int* __restrict__ nonbond_a2,
                                                      const fp_type* __restrict__ nonbond_cA,
                                                      const fp_type* __restrict__ nonbond_cB,
                                                      const int* __restrict__ nonbond_xB,
                                                      const fp_type* __restrict__ grid_maps,
                                                      const int map_index_x,
                                                      const int map_index_xy,
                                                      const int map_index_xyz) {
    const int num_atoms                = atom_counts[ligand];
    const int num_rotamers             = rotamer_counts[ligand];
    const fp_type* __restrict__ x      = x_scratch;
    const fp_type* __restrict__ y      = y_scratch;
    const fp_type* __restrict__ z      = z_scratch;
    const int atom_offset              = ligand * atoms;
    const fp_type* __restrict__ vol    = vols + atom_offset;
    const fp_type* __restrict__ solpar = solpars + atom_offset;
    const fp_type* __restrict__ charge = charges + atom_offset;
    const int* __restrict__ atom_maps  = map_offsets + atom_offset;
    const int map_size_y               = map_index_xy / map_index_x;
    const int map_size_z               = map_index_xyz / map_index_xy;
    const int last_x                   = map_index_x - 2;
    const int last_y                   = map_size_y - 2;
    const int last_z                   = map_size_z - 2;
    const fp_type* __restrict__ electro =
        grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
    const fp_type* __restrict__ desolv =
        grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

    fp_type elec_grid = 0, vdw_grid = 0, desolv_grid = 0;
    MUDOCK_PRAGMA_TEMPLATE_UNROLL(MAX_ATOMS, BLOCK_SIZE)
    for (int atom = lane; atom < MAX_ATOMS; atom += BLOCK_SIZE) {
      if (atom >= num_atoms)
        continue;
      fp_type cx = x[atom], cy = y[atom], cz = z[atom];
      if (!isfinite(cx) || !isfinite(cy) || !isfinite(cz)) {
        elec_grid += EINTCLAMP;
        vdw_grid += EINTCLAMP;
      } else if (cx < map_min[0] || cx > map_max[0] || cy < map_min[1] || cy > map_max[1] ||
                 cz < map_min[2] || cz > map_max[2]) {
        const fp_type dx      = cx - map_center[0];
        const fp_type dy      = cy - map_center[1];
        const fp_type dz      = cz - map_center[2];
        const fp_type penalty = (dx * dx + dy * dy + dz * dz) * ENERGYPENALTY;
        elec_grid += penalty;
        vdw_grid += penalty;
      } else {
        cx                = (cx - map_min[0]) * inv_spacing;
        cy                = (cy - map_min[1]) * inv_spacing;
        cz                = (cz - map_min[2]) * inv_spacing;
        const int u       = max(0, min(static_cast<int>(cx), last_x));
        const int v       = max(0, min(static_cast<int>(cy), last_y));
        const int w       = max(0, min(static_cast<int>(cz), last_z));
        const fp_type pu0 = cx - static_cast<fp_type>(u), pv0 = cy - static_cast<fp_type>(v),
                      pw0       = cz - static_cast<fp_type>(w);
        const fp_type pu[2]     = {fp_type{1} - pu0, pu0};
        const fp_type pv[2]     = {fp_type{1} - pv0, pv0};
        const fp_type pw[2]     = {fp_type{1} - pw0, pw0};
        const fp_type coeffs[8] = {pu[0] * pv[0] * pw[0],
                                   pu[0] * pv[0] * pw[1],
                                   pu[0] * pv[1] * pw[0],
                                   pu[0] * pv[1] * pw[1],
                                   pu[1] * pv[0] * pw[0],
                                   pu[1] * pv[0] * pw[1],
                                   pu[1] * pv[1] * pw[0],
                                   pu[1] * pv[1] * pw[1]};
        const int base          = FUSED_FLATTENED_3D(u, v, w, map_index_x, map_index_xy);
        const fp_type q         = charge[atom];
        elec_grid += fused_trilinear(electro + base, coeffs, map_index_x, map_index_xy) * q;
        desolv_grid += fused_trilinear(desolv + base, coeffs, map_index_x, map_index_xy) * fabsf(q);
        vdw_grid += fused_trilinear(grid_maps + atom_maps[atom] + base, coeffs, map_index_x, map_index_xy);
      }
    }

    // Previous non-unrolled intramolecular pair loop, disabled to isolate
    // its runtime cost from the grid-based score:
    // constexpr fp_type A        = mehler_solmajer::A;
    // constexpr fp_type B        = mehler_solmajer::B;
    // constexpr fp_type rk       = mehler_solmajer::rk;
    // constexpr fp_type lambda_B = mehler_solmajer::lambda_B;
    // const int nb_begin         = num_nonbond[ligand];
    // const int nb_count         = num_nonbond[ligand + 1] - nb_begin;
    // const int* nb_a1           = nonbond_a1 + nb_begin;
    // const int* nb_a2           = nonbond_a2 + nb_begin;
    // const fp_type* nb_cA       = nonbond_cA + nb_begin;
    // const fp_type* nb_cB       = nonbond_cB + nb_begin;
    // const int* nb_xB           = nonbond_xB + nb_begin;
    // fp_type elec_internal = 0, vdw_internal = 0, desolv_internal = 0;
    // if (num_rotamers > 0) {
    //   for (int i = lane; i < nb_count; i += BLOCK_SIZE) {
    //     const int a1 = nb_a1[i], a2 = nb_a2[i];
    //     const fp_type q1 = charge[a1], q2 = charge[a2];
    //     const fp_type dx = x[a1] - x[a2], dy = y[a1] - y[a2], dz = z[a1] - z[a2];
    //     const fp_type distance2  = dx * dx + dy * dy + dz * dz;
    //     const fp_type clamped    = distance2 > RMIN_ELEC_SQUARE ? distance2 : RMIN_ELEC_SQUARE;
    //     const fp_type distance   = sqrtf(clamped);
    //     const fp_type dielectric = A + B / (fp_type{1} + rk * expf(lambda_B * distance));
    //     elec_internal += q1 * q2 * ELECSCALE * autodock_parameters::coeff_estat / (distance * dielectric);
    //     const fp_type nb_desolv =
    //         vol[a2] * (solpar[a1] + qsolpar * fabsf(q1)) + vol[a1] * (solpar[a2] + qsolpar * fabsf(q2));
    //     desolv_internal +=
    //         autodock_parameters::coeff_desolv * expf(fp_type{-0.5} / sigma_square * clamped) * nb_desolv;
    //     if (clamped < nbc2 && xA_default != nb_xB[i]) {
    //       const fp_type log_distance = logf(distance);
    //       const fp_type rA           = expf(static_cast<fp_type>(xA_default) * log_distance);
    //       const fp_type rB           = expf(static_cast<fp_type>(nb_xB[i]) * log_distance);
    //       fp_type e                  = nb_cA[i] / rA - nb_cB[i] / rB;
    //       e                          = EINTCLAMP < e ? EINTCLAMP : e;
    //       vdw_internal += e;
    //     }
    //   }
    // }
    // fp_type elec_energy = 0, vdw_energy = 0, desolv_energy = 0;
    // elec_energy    = elec_internal;
    // vdw_energy     = vdw_internal;
    // desolv_energy  = desolv_internal;
    fp_type energy = vdw_grid + elec_grid + desolv_grid;
    for (int offset = BLOCK_SIZE / 2; offset > 0; offset /= 2)
      energy += __shfl_down_sync(0xffffffff, energy, offset);
    if (lane == 0)
      energy += static_cast<fp_type>(num_rotamers) * autodock_parameters::coeff_tors;
    return __shfl_sync(0xffffffff, energy, 0);
  }

  template<int MAX_ATOMS, int MAX_ROTAMERS>
  __global__ void FUSED_GENETIC_LAUNCH_BOUNDS(MAX_ATOMS)
      genetic_fused_kernel(const int ligands,
                           const int atoms,
                           const int population,
                           const int generations,
                           const int tournament_length,
                           const fp_type mutation_probability,
                           const std::size_t seed,
                           fp_type* __restrict__ current_data,
                           fp_type* __restrict__ next_data,
                           fp_type* __restrict__ scores_data,
                           fp_type* __restrict__ best_scores,
                           fp_type* __restrict__ best_chromosomes,
                           const int* __restrict__ atom_counts,
                           const int* __restrict__ rotamer_counts,
                           const fp_type* __restrict__ x_coords,
                           const fp_type* __restrict__ y_coords,
                           const fp_type* __restrict__ z_coords,
                           const fp_type* __restrict__ map_min,
                           const fp_type* __restrict__ map_max,
                           const fp_type* __restrict__ map_center,
                           const int* __restrict__ fragments,
                           const int* __restrict__ fragments_start,
                           const int* __restrict__ frag_indices_start,
                           const int* __restrict__ frag_start_data,
                           const int* __restrict__ frag_stop_data,
                           const int* __restrict__ num_nonbond,
                           const fp_type* __restrict__ vols,
                           const fp_type* __restrict__ solpars,
                           const fp_type* __restrict__ charges,
                           const int* __restrict__ map_offsets,
                           const int* __restrict__ nonbond_a1,
                           const int* __restrict__ nonbond_a2,
                           const fp_type* __restrict__ nonbond_cA,
                           const fp_type* __restrict__ nonbond_cB,
                           const int* __restrict__ nonbond_xB,
                           const fp_type* __restrict__ grid_maps,
                           const int map_index_x,
                           const int map_index_xy,
                           const int map_index_xyz) {
    const int ligand = blockIdx.x;
    if (ligand >= ligands)
      return;
    const int lane                   = threadIdx.x % BLOCK_SIZE;
    const int warp                   = threadIdx.x / BLOCK_SIZE;
    const int thread                 = threadIdx.x;
    const int n_atoms                = atom_counts[ligand];
    const int n_rotamers             = rotamer_counts[ligand];
    constexpr int gene_stride = 6 + MAX_ROTAMERS;
    const int ligand_gene_offset = ligand * population * gene_stride;
    const int ligand_score_offset = ligand * population;
    fp_type* __restrict__ current = current_data + ligand_gene_offset;
    fp_type* __restrict__ next    = next_data + ligand_gene_offset;
    fp_type* __restrict__ scores   = scores_data + ligand_score_offset;
    __shared__ fp_type pose_x[FUSED_WARPS_PER_BLOCK][MAX_ATOMS];
    __shared__ fp_type pose_y[FUSED_WARPS_PER_BLOCK][MAX_ATOMS];
    __shared__ fp_type pose_z[FUSED_WARPS_PER_BLOCK][MAX_ATOMS];

    curandState rng;
    const unsigned long long tid = static_cast<unsigned long long>(ligand) * FUSED_BLOCK_SIZE + thread;
    curand_init(seed + tid, tid, 0, &rng);

    for (int i = thread; i < population; i += FUSED_BLOCK_SIZE) {
      fp_type* genes = current + i * gene_stride;
      scores[i] = std::numeric_limits<fp_type>::infinity();
      for (int gene = 0; gene < 3; ++gene) genes[gene] = fused_random(rng, -45, 45) * fp_type{0.2};
      for (int gene = 3; gene < 6 + n_rotamers; ++gene)
        genes[gene] = fused_random(rng, -45, 45) * fp_type{4};
    }
    __syncthreads();

    const int atom_base                      = ligand * atoms;
    const fp_type* __restrict__ lx           = x_coords + atom_base;
    const fp_type* __restrict__ ly           = y_coords + atom_base;
    const fp_type* __restrict__ lz           = z_coords + atom_base;
    const int* __restrict__ ligand_fragments = fragments + fragments_start[ligand];
    const int* __restrict__ frag_start       = frag_start_data + frag_indices_start[ligand];
    const int* __restrict__ frag_stop        = frag_stop_data + frag_indices_start[ligand];

    for (int generation = 0; generation < generations; ++generation) {
      for (int pose_base = 0; pose_base < population; pose_base += FUSED_WARPS_PER_BLOCK) {
        const int pose = pose_base + warp;
        if (pose < population) {
          fp_type* __restrict__ x = pose_x[warp];
          fp_type* __restrict__ y = pose_y[warp];
          fp_type* __restrict__ z = pose_z[warp];
          MUDOCK_PRAGMA_TEMPLATE_UNROLL(MAX_ATOMS, BLOCK_SIZE)
          for (int atom = lane; atom < MAX_ATOMS; atom += BLOCK_SIZE) {
            if (atom < n_atoms) {
              x[atom] = lx[atom];
              y[atom] = ly[atom];
              z[atom] = lz[atom];
            }
          }
          __syncwarp();
          translate_molecule_cuda<MAX_ATOMS>(x,
                                             y,
                                             z,
                                             current[pose * gene_stride],
                                             current[pose * gene_stride + 1],
                                             current[pose * gene_stride + 2],
                                             n_atoms);
          __syncwarp();
          rotate_molecule_cuda<MAX_ATOMS>(x,
                                          y,
                                          z,
                                          current[pose * gene_stride + 3],
                                          current[pose * gene_stride + 4],
                                          current[pose * gene_stride + 5],
                                          n_atoms);
          __syncwarp();
          MUDOCK_PRAGMA_TEMPLATE_UNROLL(MAX_ROTAMERS / FUSED_ROTAMER_UNROLL_DIVISOR, 1)
          for (int rot = 0; rot < MAX_ROTAMERS; ++rot) {
            if (rot < n_rotamers) {
              rotate_fragment_cuda<MAX_ATOMS>(x,
                                              y,
                                              z,
                                              ligand_fragments + rot * n_atoms,
                                              frag_start[rot],
                                              frag_stop[rot],
                                              current[pose * gene_stride + 6 + rot],
                                              n_atoms);
              __syncwarp();
            }
          }
          const fp_type score = fused_score_pose<MAX_ATOMS>(ligand,
                                                            lane,
                                                            atoms,
                                                            atom_counts,
                                                            rotamer_counts,
                                                            x,
                                                            y,
                                                            z,
                                                            map_min,
                                                            map_max,
                                                            map_center,
                                                            vols,
                                                            solpars,
                                                            charges,
                                                            map_offsets,
                                                            num_nonbond,
                                                            nonbond_a1,
                                                            nonbond_a2,
                                                            nonbond_cA,
                                                            nonbond_cB,
                                                            nonbond_xB,
                                                            grid_maps,
                                                            map_index_x,
                                                            map_index_xy,
                                                            map_index_xyz);
          if (lane == 0)
            scores[pose] = score;
        }
      }
      __syncthreads();

      if (generation + 1 < generations) {
        for (int child = thread; child < population; child += FUSED_BLOCK_SIZE) {
          const int p1    = fused_select(rng, population, tournament_length, scores);
          const int p2    = fused_select(rng, population, tournament_length, scores);
          const int split = fused_random_int(rng, 0, 6 + n_rotamers);
          fp_type* child_genes = next + child * gene_stride;
          const fp_type* parent1_genes = current + p1 * gene_stride;
          const fp_type* parent2_genes = current + p2 * gene_stride;
          for (int gene = 0; gene < 6 + n_rotamers; ++gene)
            child_genes[gene] = gene < split ? parent1_genes[gene] : parent2_genes[gene];
          for (int gene = 0; gene < 3; ++gene)
            if (fused_random(rng, 0, 1) < mutation_probability)
              child_genes[gene] += fused_random(rng, -10, 10) * fp_type{0.2};
          for (int gene = 3; gene < 6 + n_rotamers; ++gene)
            if (fused_random(rng, 0, 1) < mutation_probability)
              child_genes[gene] += fused_random(rng, -10, 10) * fp_type{4};
        }
        __syncthreads();
        fp_type* swap = current;
        current          = next;
        next             = swap;
      }
    }

    if (warp == 0) {
      int best_index     = lane;
      fp_type best_score = lane < population ? scores[lane] : std::numeric_limits<fp_type>::infinity();
      for (int i = lane + BLOCK_SIZE; i < population; i += BLOCK_SIZE) {
        if (scores[i] < best_score) {
          best_score = scores[i];
          best_index = i;
        }
      }
      for (int delta = BLOCK_SIZE / 2; delta > 0; delta /= 2) {
        const fp_type other_score = __shfl_down_sync(0xffffffff, best_score, delta);
        const int other_index     = __shfl_down_sync(0xffffffff, best_index, delta);
        if (other_score < best_score) {
          best_score = other_score;
          best_index = other_index;
        }
      }
      if (lane == 0) {
        best_scores[ligand] = best_score;
        const fp_type* best_genes = current + best_index * gene_stride;
        fp_type* output_genes = best_chromosomes + ligand * gene_stride;
        for (int gene = 0; gene < 6 + n_rotamers; ++gene)
          output_genes[gene] = best_genes[gene];
      }
    }
  }

  void genetic_fused_launcher<queue_cuda>::operator()() {
    void* args[] = {&data.ligands,
                    &data.atoms,
                    &data.population,
                    &data.generations,
                    &data.tournament_length,
                    &data.mutation_probability,
                    &data.seed,
                    &data.current,
                    &data.next,
                    &data.scores,
                    &data.best_scores,
                    &data.best_chromosomes,
                    &data.num_atoms,
                    &data.num_rotamers,
                    &data.x_coords,
                    &data.y_coords,
                    &data.z_coords,
                    &data.map_min,
                    &data.map_max,
                    &data.map_center,
                    &data.fragments,
                    &data.fragments_start,
                    &data.frag_indices_start,
                    &data.frag_start,
                    &data.frag_stop,
                    &data.num_nonbond,
                    &data.vols,
                    &data.solpars,
                    &data.charges,
                    &data.map_offsets,
                    &data.nonbond_a1,
                    &data.nonbond_a2,
                    &data.nonbond_cA,
                    &data.nonbond_cB,
                    &data.nonbond_xB,
                    &data.grid_maps,
                    &data.map_index_x,
                    &data.map_index_xy,
                    &data.map_index_xyz};
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_cluster) {
          constexpr int max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_cluster];
          constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_rotamer_clusters(), 1>(
              [&](const auto rotamer_cluster) {
                constexpr int max_rotamers =
                    reorder_buffer<static_molecule>::rotamers_clusters[rotamer_cluster];
                q->launch_kernel((void*) genetic_fused_kernel<max_atoms, max_rotamers>,
                                 args,
                                 data.ligands,
                                 FUSED_BLOCK_SIZE);
              },
              data.max_rotamers,
              reorder_buffer<static_molecule>::rotamers_clusters.data());
        },
        data.atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }

  batch_multiple get_genetic_fused_batch_multiple_cuda(const int device_id,
                                                        const int atoms,
                                                        const int max_rotamers) {
    batch_multiple multiple{};
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_cluster) {
          constexpr int max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_cluster];
          constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_rotamer_clusters(), 1>(
              [&](const auto rotamer_cluster) {
                constexpr int max_rotamers_bucket =
                    reorder_buffer<static_molecule>::rotamers_clusters[rotamer_cluster];
                multiple = get_kernel_batch_multiple_cuda<
                    genetic_fused_kernel<max_atoms, max_rotamers_bucket>>(device_id, FUSED_BLOCK_SIZE);
              },
              max_rotamers,
              reorder_buffer<static_molecule>::rotamers_clusters.data());
        },
        atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());
    return multiple;
  }
} // namespace mudock
